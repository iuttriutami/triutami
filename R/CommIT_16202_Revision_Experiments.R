# =============================================================================
# CommIT 16202 - REVISION ROUND 1: CONSOLIDATED & REPRODUCIBLE EXPERIMENTS
# "Metaheuristic-Optimized K-Medoids Clustering of Indonesian Provinces Based on
#  Environmental Quality Index Components" (EQI / IKLH 2023, 38 provinces)
# -----------------------------------------------------------------------------
# This single script replaces KMedoids_metaheuristic_COMPLETE_EN.R and
# CommIT_EqualBudget_Original_R.R for the revised manuscript (v10).
#
# What it produces (all in DIR_OUT):
#   E0  preprocessing: Mahalanobis D^2 (cutoff chi2_{p;0.999}), VIF, correlations
#   E1  PAM baseline (cluster::pam, BUILD+SWAP) for k = 2..10           [R3]
#   E2  exhaustive enumeration (global optimum) for k = 2..6            [R2]
#   E3  EQUAL-BUDGET benchmark: GA, PSO, CSA, GWO + random search (RS),
#       exactly 3,000 fitness evaluations per run, 30 runs, k = 2..10   [R3, R1]
#       -> Table 2, Table 5, Friedman + Wilcoxon-Holm, convergence plot
#   E4  NOMINAL-BUDGET runs (original settings) with EXACT evaluation counts [R3]
#   E5  sensitivity to population size under the equal budget              [R2]
#   E6  final k = 3 (and k = 2) partitions, per-object silhouettes, profile [R2, R3]
#   E7  DBI / CH across k                                                  [R2]
#   E8  leave-one-province-out stability (ARI)                             [R2]
#   E9  Silhouette-vs-k figure                                             [R2]
#   E10 corrected choropleth map with strict name matching + Java inset   [R1, R2]
#   placeholders.csv  -> every [[TOKEN]] used in the revised manuscript
#
# Two corrections relative to the v9 scripts (please mention in the response):
#   (a) the medoid set reported for each run is the set that was ACTUALLY
#       evaluated when the best fitness was obtained; v9 re-decoded the final
#       continuous solution, and the stochastic duplicate repair could then
#       return a different (worse) medoid set (e.g., PSO k = 5 SC = 0.0560).
#   (b) fitness evaluations are COUNTED, not inferred from population x
#       iterations (GA::ga does not re-evaluate unchanged elites).
#
# Fitness: mean Silhouette with nearest-medoid assignment. A vectorised
# implementation is used for speed and is verified to equal
# cluster::silhouette() to 1e-10 before any experiment is run.
# =============================================================================

# ---- 0. SETTINGS --------------------------------------------------------------
SEED          <- 42L
N_RUN         <- 30L            # independent runs; run r uses seed SEED + r (43..72)
BUDGET        <- 3000L          # exact fitness evaluations per run (equal budget)
K_FINAL       <- 3L
K_RANGE       <- 2:10
EXHAUSTIVE_K  <- 2:6            # C(38,6) = 2,760,681 sets; ~2-6 min with fast fitness
META          <- c("GA", "PSO", "CSA", "GWO")
ALGS_EB       <- c(META, "RS")  # RS = pure random search, same budget (reference)
RUN_SENSITIVITY <- TRUE
RUN_LOO         <- TRUE

ORIGINAL_NAME <- "Komponen Penyusun Indeks Kualitas Lingkungan Hidup Menurut Provinsi, 2023.csv"
PATH_CSV <- Sys.getenv("EQI_CSV_PATH", "")   # or set directly, e.g. "C:/.../file.csv"
PATH_SHP <- Sys.getenv("EQI_SHP_PATH", "")   # e.g. "indo_by_prov_2023.shp"; "" = search
DIR_OUT  <- if (dir.exists("/kaggle/working")) "/kaggle/working/commIT_revision" else "commIT_revision"
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

pkgs <- c("cluster", "GA", "pso", "clusterCrit", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, repos = "https://cloud.r-project.org")
HAS_SF <- requireNamespace("sf", quietly = TRUE)
suppressPackageStartupMessages({ library(cluster); library(ggplot2) })

T0_ALL <- proc.time()[["elapsed"]]
TOK <- new.env()                                   # manuscript placeholders
tok <- function(key, value) assign(key, as.character(value), envir = TOK)
f4  <- function(x) ifelse(is.na(x), "\u2013", sprintf("%.4f", x))
f3  <- function(x) sprintf("%.3f", x)
f2  <- function(x) sprintf("%.2f", x)
fp  <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)))
fpct <- function(x) ifelse(is.na(x), "\u2013", sprintf("%.0f%%", 100 * x))
fint <- function(x) format(round(x), big.mark = ",", trim = TRUE)
say <- function(...) cat(sprintf(...), "\n")
wcsv <- function(x, f) write.csv(x, file.path(DIR_OUT, f), row.names = FALSE)

# ---- 1. DATA -----------------------------------------------------------------
find_file <- function(explicit, exact_name, pattern) {
  if (nzchar(explicit)) { if (!file.exists(explicit)) stop("Not found: ", explicit); return(explicit) }
  for (root in c("/kaggle/input", "/mnt/data", getwd())) {
    if (!dir.exists(root)) next
    cand <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (is.null(exact_name)) { if (length(cand)) return(cand[1]) else next }
    hit <- cand[basename(cand) == exact_name]
    if (!length(hit)) hit <- cand[grepl("Kualitas Lingkungan", basename(cand), ignore.case = TRUE)]
    if (!length(hit)) hit <- Filter(function(f) {          # any CSV whose header has the EQI columns
      h <- tryCatch(readLines(f, n = 1, warn = FALSE), error = function(e) "")
      length(h) && grepl("Indeks Kualitas Udara", h, fixed = TRUE) }, cand)
    if (length(hit)) return(hit[1])
  }
  NULL
}
find_shp <- function(explicit) {   # never pick a regional (e.g. Jawa Tengah) shapefile by accident
  if (nzchar(explicit)) return(if (file.exists(explicit)) explicit else NULL)
  for (root in c("/kaggle/input", "/mnt/data", getwd())) {
    if (!dir.exists(root)) next
    cand <- list.files(root, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
    cand <- cand[grepl("prov|indo", cand, ignore.case = TRUE) & !grepl("jateng|kab|kota|desa|kec", cand, ignore.case = TRUE)]
    if (length(cand)) return(cand[1])
  }
  NULL
}
path <- find_file(PATH_CSV, ORIGINAL_NAME, "\\.csv$")
if (is.null(path)) stop("EQI CSV not found. Set EQI_CSV_PATH or PATH_CSV.")
df <- read.csv(path, header = TRUE, check.names = FALSE, na.strings = c("", "NA", "-"))
names(df)[1] <- "Provinsi"
df <- df[!toupper(trimws(df$Provinsi)) %in% c("INDONESIA", ""), , drop = FALSE]
req <- c("Indeks Kualitas Udara", "Indeks Kualitas Air",
         "Indeks Kualitas Tutupan Lahan (IKTL)", "Indeks Kualitas Lahan (IKL)")
stopifnot(all(req %in% names(df)))
provinsi <- trimws(as.character(df$Provinsi))
data_num <- df[, req]
for (j in seq_along(req)) data_num[[j]] <- as.numeric(data_num[[j]])
SHORT <- c("AQI", "WQI", "LCQI", "LQI")
names(data_num) <- SHORT
stopifnot(nrow(data_num) == 38L, !anyNA(data_num), !anyDuplicated(provinsi))
say("Data: %s | %d provinces x %d indicators", basename(path), nrow(data_num), ncol(data_num))
n <- nrow(data_num); p <- ncol(data_num)
PU <- toupper(gsub("\\s+", " ", provinsi))

# ---- E0. PREPROCESSING DIAGNOSTICS ------------------------------------------
cat("\n=== E0 Missing values per column ===\n"); print(colSums(is.na(data_num)))
D2 <- mahalanobis(as.matrix(data_num), colMeans(data_num), cov(data_num))   # squared distance
CUT <- qchisq(0.999, df = p)
mah <- data.frame(Province = provinsi, D2 = D2, d = sqrt(D2), Outlier = D2 > CUT)
mah <- mah[order(-mah$D2), ]; mah$Rank <- seq_len(n)
wcsv(transform(mah, D2 = round(D2, 4), d = round(d, 4)), "E0_mahalanobis.csv")
say("Mahalanobis: cutoff chi2(%d; 0.999) = %.2f ; max D2 = %.4f (%s) ; flagged = %d",
    p, CUT, mah$D2[1], mah$Province[1], sum(mah$Outlier))
tok("MAHAL_CUT", f2(CUT)); tok("MAHAL_MAX_D2", f4(mah$D2[1])); tok("MAHAL_MAX_PROV", mah$Province[1])
tok("MAHAL_NOUT", sum(mah$Outlier))
tok("MAHAL_OUT_LIST", if (any(mah$Outlier)) paste(mah$Province[mah$Outlier], collapse = ", ") else "none")
tok("MAHAL_TOP3", paste(sprintf("%s (D\u00b2 = %s)", mah$Province[1:3], f2(mah$D2[1:3])), collapse = ", "))
tok("MAHAL_RANK_JKT", mah$Rank[toupper(mah$Province) == "DKI JAKARTA"])
tok("MAHAL_RANK_BTN", mah$Rank[toupper(mah$Province) == "BANTEN"])
if (mah$D2[1] < 6) message("NOTE: max D2 is small; v9 reported 4.5341 - check whether v9 reported d or D2.")

vif <- sapply(seq_len(p), function(j)
  1 / (1 - summary(lm(data_num[[j]] ~ ., data = data_num[, -j, drop = FALSE]))$r.squared))
wcsv(data.frame(Variable = SHORT, VIF = round(vif, 4)), "E0_vif.csv")
wcsv(round(cor(data_num), 3), "E0_correlation.csv")
tok("VIF_MAX", f4(max(vif))); tok("COR_MAX", f3(max(abs(cor(data_num)[upper.tri(diag(p))]))))

X    <- scale(data_num)            # z-scores (sample SD)
DIST <- dist(X)
DMAT <- as.matrix(DIST)

# ---- HELPERS -----------------------------------------------------------------
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

# verify fast fitness == cluster::silhouette
set.seed(1)
chk <- unlist(lapply(K_RANGE, function(k) replicate(300, {
  med <- sample(n, k); lab <- labels_of(med)
  abs(fast_sil_lab(lab) - mean(cluster::silhouette(lab, DIST)[, 3]))
})))
stopifnot(max(chk) < 1e-10)
say("Fast Silhouette verified against cluster::silhouette (max |diff| = %.1e)", max(chk))
ref_med <- match(c("RIAU", "DKI JAKARTA", "BANTEN"), PU)
REF_SC <- sil_med(ref_med)
say("v9 reference partition (Riau | DKI Jakarta | Banten): SC = %.6f", REF_SC)
if (abs(REF_SC - 0.4629) > 5e-5) warning("Reference SC differs from v9 (0.4629): check data file!")

# ---- E1. PAM BASELINE ------------------------------------------------------------
pam_tab <- do.call(rbind, lapply(K_RANGE, function(k) {
  pm <- cluster::pam(DIST, k = k, diss = TRUE)     # BUILD initialisation + SWAP, deterministic
  data.frame(k = k, PAM_SC = sil_med(pm$id.med), PAM_obj = pm$objective[["swap"]] * n,
             PAM_sizes = paste(sort(tabulate(pm$clustering)), collapse = "-"))
}))
wcsv(pam_tab, "E1_pam_baseline.csv"); print(pam_tab)
for (k in K_RANGE) tok(sprintf("PAM_K%d", k), f4(pam_tab$PAM_SC[pam_tab$k == k]))

# ---- E2. EXHAUSTIVE ENUMERATION -------------------------------------------------
glob <- list()
for (k in EXHAUSTIVE_K) {
  t0 <- proc.time()[["elapsed"]]
  cmb <- utils::combn(n, k)
  sc <- vapply(seq_len(ncol(cmb)), function(j) sil_med(cmb[, j]), numeric(1))
  b <- which.max(sc)
  glob[[as.character(k)]] <- list(k = k, n_sets = ncol(cmb), best = sc[b], med = cmb[, b],
    n_opt = sum(sc > sc[b] - 1e-12), secs = proc.time()[["elapsed"]] - t0,
    pam_pct = mean(sc <= pam_tab$PAM_SC[pam_tab$k == k] + 1e-12))
  say("Exhaustive k=%d: %s sets, global max SC = %.6f (%s) [%.0f s]", k, fint(ncol(cmb)),
      sc[b], paste(provinsi[cmb[, b]], collapse = " | "), glob[[as.character(k)]]$secs)
  rm(cmb, sc); invisible(gc())
}
glob_tab <- do.call(rbind, lapply(glob, function(g) data.frame(k = g$k, sets = g$n_sets,
  global_SC = g$best, n_optimal_sets = g$n_opt, medoids = paste(provinsi[g$med], collapse = " | "),
  PAM_percentile = g$pam_pct, seconds = round(g$secs, 1))))
wcsv(glob_tab, "E2_exhaustive_global_optimum.csv")
GLOBAL <- setNames(glob_tab$global_SC, glob_tab$k)
for (k in EXHAUSTIVE_K) { tok(sprintf("GLOBAL_K%d", k), f4(GLOBAL[[as.character(k)]]))
  tok(sprintf("NSETS_K%d", k), fint(glob_tab$sets[glob_tab$k == k])) }
tok("EXH_KMAX", max(EXHAUSTIVE_K))
tok("EXH_SECS_KMAX", sprintf("%.0f", glob_tab$seconds[glob_tab$k == max(EXHAUSTIVE_K)]))
if (abs(GLOBAL[["3"]] - REF_SC) > 1e-9) warning("Global k=3 optimum differs from the v9 reference partition!")

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

# ---- E3. EQUAL-BUDGET BENCHMARK ------------------------------------------------------
say("\n=== E3 Equal-budget benchmark: %d evals x %d runs x %d algorithms x k = %d..%d ===",
    BUDGET, N_RUN, length(ALGS_EB), min(K_RANGE), max(K_RANGE))
rows <- list(); traces3 <- list(); best_part <- list()
for (k in K_RANGE) {
  tgt <- if (as.character(k) %in% names(GLOBAL)) GLOBAL[[as.character(k)]] else NA_real_
  for (r in seq_len(N_RUN)) for (alg in ALGS_EB) {
    o <- run_one(alg, k, SEED + r, target = tgt)
    rows[[length(rows) + 1]] <- data.frame(k = k, run = r, seed = SEED + r, algorithm = alg,
      best_SC = o$best, evals = o$calls, evals_to_opt = o$hit, success = !is.na(o$hit),
      secs_to_budget = o$secs, medoids = paste(provinsi[o$med], collapse = " | "),
      sizes = paste(sort(tabulate(labels_of(o$med), k)), collapse = "-"))
    if (k == K_FINAL) traces3[[length(traces3) + 1]] <- data.frame(algorithm = alg, run = r,
      evaluation = seq_len(BUDGET), best_SC = o$trace)
    key <- as.character(k)
    if (is.null(best_part[[key]]) || o$best > best_part[[key]]$best + 1e-12) best_part[[key]] <- o
  }
  say("k = %d done", k)
}
eb <- do.call(rbind, rows); wcsv(eb, "E3_equal_budget_runs.csv")
if (any(is.na(match(names(GLOBAL), unique(eb$k))))) stop("k mismatch")
eb_sum <- do.call(rbind, lapply(split(eb, list(eb$algorithm, eb$k), drop = TRUE), function(d)
  data.frame(k = d$k[1], algorithm = d$algorithm[1], mean_SC = mean(d$best_SC), sd_SC = sd(d$best_SC),
    median_SC = median(d$best_SC), min_SC = min(d$best_SC), max_SC = max(d$best_SC),
    success_rate = if (all(is.na(d$evals_to_opt)) && !(as.character(d$k[1]) %in% names(GLOBAL))) NA else mean(d$success),
    median_evals_to_opt = if (any(d$success)) median(d$evals_to_opt[d$success]) else NA,
    mean_secs = mean(d$secs_to_budget), sd_secs = sd(d$secs_to_budget))))
eb_sum <- eb_sum[order(eb_sum$k, match(eb_sum$algorithm, ALGS_EB)), ]
wcsv(eb_sum, "E3_equal_budget_summary.csv"); print(eb_sum, row.names = FALSE)
G <- function(k, a, col) eb_sum[eb_sum$k == k & eb_sum$algorithm == a, col]

# Table 2 (manuscript): PAM | global optimum | mean best SC over 30 equal-budget runs
t2 <- data.frame(k = K_RANGE, PAM = pam_tab$PAM_SC,
  Global = sapply(K_RANGE, function(k) if (as.character(k) %in% names(GLOBAL)) GLOBAL[[as.character(k)]] else NA))
for (a in ALGS_EB) t2[[a]] <- sapply(K_RANGE, function(k) G(k, a, "mean_SC"))
wcsv(cbind(k = t2$k, as.data.frame(lapply(t2[-1], f4))), "TABLE2_for_manuscript.csv")
for (k in K_RANGE) {
  tok(sprintf("T2_PAM_K%d", k), f4(t2$PAM[t2$k == k])); tok(sprintf("T2_GLOB_K%d", k), f4(t2$Global[t2$k == k]))
  for (a in ALGS_EB) tok(sprintf("T2_%s_K%d", a, k), f4(t2[[a]][t2$k == k]))
}
# where does PAM beat the mean of a metaheuristic?
beat <- do.call(rbind, lapply(K_RANGE, function(k) { m <- sapply(META, function(a) G(k, a, "mean_SC"))
  if (any(m < pam_tab$PAM_SC[pam_tab$k == k])) data.frame(k = k, algs = paste(META[m < pam_tab$PAM_SC[pam_tab$k == k]], collapse = "/")) }))
tok("PAM_BEATS_SENTENCE", if (is.null(beat)) "no metaheuristic mean fell below the PAM baseline at any k" else
  paste0("the PAM baseline exceeded the mean equal-budget Silhouette of ",
         paste(sprintf("%s at k = %d", beat$algs, beat$k), collapse = "; ")))
small_k <- intersect(2:4, K_RANGE)
tok("GAIN_K3", f4(GLOBAL[["3"]] - pam_tab$PAM_SC[pam_tab$k == 3]))

# Table 5 (equal budget, k = 3) + stats
for (a in ALGS_EB) {
  tok(sprintf("EB3_MEAN_%s", a), f4(G(K_FINAL, a, "mean_SC"))); tok(sprintf("EB3_SD_%s", a), f4(G(K_FINAL, a, "sd_SC")))
  tok(sprintf("EB3_SR_%s", a), fpct(G(K_FINAL, a, "success_rate")))
  me <- G(K_FINAL, a, "median_evals_to_opt"); tok(sprintf("EB3_E2O_%s", a), ifelse(is.na(me), "\u2013", fint(me)))
  tok(sprintf("EB3_SEC_%s", a), f2(G(K_FINAL, a, "mean_secs")))
}
wide_of <- function(k, algs) sapply(algs, function(a) { d <- eb[eb$k == k & eb$algorithm == a, ]; round(d$best_SC[order(d$run)], 10) })  # drop float noise between equivalent medoid sets
stat_k <- function(k, algs = META) {
  W <- wide_of(k, algs)
  fr <- if (all(abs(W - W[1, 1]) < 1e-12)) list(statistic = NA, p.value = NA) else
    tryCatch(friedman.test(W), error = function(e) list(statistic = NA, p.value = NA))
  pr <- combn(algs, 2, simplify = FALSE)
  wx <- do.call(rbind, lapply(pr, function(ab) { d <- W[, ab[1]] - W[, ab[2]]
    pv <- if (all(abs(d) < 1e-12)) 1 else suppressWarnings(wilcox.test(W[, ab[1]], W[, ab[2]],
      paired = TRUE, exact = FALSE, correct = FALSE)$p.value)
    data.frame(k = k, A = ab[1], B = ab[2], mean_diff = mean(d), nonzero = sum(abs(d) > 1e-12), p_raw = pv) }))
  wx$p_holm <- p.adjust(wx$p_raw, "holm")
  list(chi = unname(fr$statistic), p = unname(fr$p.value), wx = wx)
}
ST <- lapply(K_RANGE, stat_k); names(ST) <- K_RANGE
wcsv(data.frame(k = K_RANGE, chi2 = sapply(ST, `[[`, "chi"), df = 3, p = sapply(ST, `[[`, "p")), "E3_friedman_by_k.csv")
wcsv(do.call(rbind, lapply(ST, `[[`, "wx")), "E3_wilcoxon_holm_by_k.csv")
S3 <- ST[[as.character(K_FINAL)]]
tok("EB3_FR_CHI", ifelse(is.na(S3$chi), "n/a", f3(S3$chi))); tok("EB3_FR_P", fp(S3$p))
tok("EB3_WX_MINP", fp(min(S3$wx$p_holm)))
sig <- S3$wx[S3$wx$p_holm < 0.05, ]
tok("EB3_STATS_SENTENCE", if (is.na(S3$p)) {
  "all four optimizers returned the same Silhouette in every run, so the rank-based tests are undefined and no difference can be claimed"
} else if (S3$p < 0.05) {
  paste0("the Friedman test detected a statistically significant difference among the four optimizers (\u03c7\u00b2 = ",
    f3(S3$chi), ", df = 3, p ", fp(S3$p), ")", if (nrow(sig)) paste0(", and the Holm-adjusted Wilcoxon signed-rank tests located it in the pairs ",
    paste(sprintf("%s\u2013%s", sig$A, sig$B), collapse = ", ")) else ", although no Holm-adjusted pairwise Wilcoxon test was significant")
} else {
  paste0("the Friedman test did not detect a statistically significant difference among the four optimizers (\u03c7\u00b2 = ",
    f3(S3$chi), ", df = 3, p ", fp(S3$p), "), and no Holm-adjusted pairwise Wilcoxon test was significant at \u03b1 = 0.05")
})
# ranking by success rate, then median evaluations-to-optimum, then runtime
rk <- eb_sum[eb_sum$k == K_FINAL & eb_sum$algorithm %in% META, ]
rk <- rk[order(-rk$success_rate, rk$median_evals_to_opt, rk$mean_secs), ]
tok("EB3_RANKING", paste(rk$algorithm, collapse = " > "))
tok("EB3_TOP", rk$algorithm[1])
tok("EB3_ALL_OPT_ALGS", { a <- rk$algorithm[rk$success_rate == 1]; if (length(a)) paste(a, collapse = ", ") else "none of the optimizers" })
if (rk$algorithm[1] != "CSA") warning("AUTO-CHECK: CSA is NOT the top-ranked optimizer under the equal budget. ",
  "Revise the sentences that retain CSA as the preferred optimizer (search for EB3_TOP in the manuscript).")
tok("EB3_RS_SR", fpct(G(K_FINAL, "RS", "success_rate")))
for (k in intersect(c(2, 4, 5, 6), K_RANGE)) for (a in META) tok(sprintf("EB_SR_%s_K%d", a, k), fpct(G(k, a, "success_rate")))
tok("EB_SR_RANGE_K6", { s <- sapply(META, function(a) G(6, a, "success_rate")); sprintf("%s\u2013%s", fpct(min(s)), fpct(max(s))) })
tok("EB_RS_SR_K6", fpct(G(6, "RS", "success_rate")))

# convergence by evaluations (k = 3)
tr <- do.call(rbind, traces3)
mtr <- aggregate(best_SC ~ algorithm + evaluation, data = tr, FUN = mean)
wcsv(mtr, "E3_mean_convergence_k3.csv")
mtr$algorithm <- factor(mtr$algorithm, levels = ALGS_EB)
pc <- ggplot(mtr, aes(evaluation, best_SC, colour = algorithm, linetype = algorithm)) +
  geom_line(linewidth = 0.8) + geom_hline(yintercept = GLOBAL[["3"]], linetype = 3, colour = "grey30") +
  scale_linetype_manual(values = c(GA = 1, PSO = 1, CSA = 1, GWO = 1, RS = 2)) +
  scale_colour_manual(values = c(GA = "#1b9e77", PSO = "#d95f02", CSA = "#7570b3", GWO = "#e7298a", RS = "grey45")) +
  labs(x = "Cumulative fitness evaluations", y = "Mean best-so-far Silhouette (30 runs)", colour = NULL, linetype = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(DIR_OUT, "Fig_convergence_equal_budget_k3.png"), pc, width = 6.5, height = 4, dpi = 300)
bx <- eb[eb$k == K_FINAL, ]; bx$algorithm <- factor(bx$algorithm, levels = ALGS_EB)
ggsave(file.path(DIR_OUT, "FigS_boxplot_equal_budget_k3.png"),
  ggplot(bx, aes(algorithm, best_SC)) + geom_boxplot() + geom_jitter(width = 0.1, alpha = 0.4) +
    geom_hline(yintercept = GLOBAL[["3"]], linetype = 3) + labs(x = NULL, y = "Best Silhouette per run") +
    theme_bw(), width = 6, height = 4, dpi = 300)

# ---- E4. NOMINAL-BUDGET RUNS (v9 settings) with exact evaluation counts --------------
nom <- do.call(rbind, lapply(seq_len(N_RUN), function(r) do.call(rbind, lapply(META, function(a) {
  o <- run_one(a, K_FINAL, SEED + r, budget = 1e7, pop = NOMINAL[[a]]$pop,
               maxiter = NOMINAL[[a]]$maxiter, target = GLOBAL[["3"]])
  data.frame(run = r, algorithm = a, best_SC = o$best, evals = o$calls, success = !is.na(o$hit), secs = o$secs)
}))))
wcsv(nom, "E4_nominal_budget_runs_k3.csv")
nom_sum <- do.call(rbind, lapply(META, function(a) { d <- nom[nom$algorithm == a, ]
  data.frame(algorithm = a, mean_evals = mean(d$evals), mean_SC = mean(d$best_SC), sd_SC = sd(d$best_SC),
             success_rate = mean(d$success), mean_secs = mean(d$secs)) }))
wcsv(nom_sum, "E4_nominal_budget_summary_k3.csv"); print(nom_sum, row.names = FALSE)
for (a in META) { d <- nom_sum[nom_sum$algorithm == a, ]
  tok(sprintf("NOM_EVALS_%s", a), fint(d$mean_evals)); tok(sprintf("NOM_MEAN_%s", a), f4(d$mean_SC))
  tok(sprintf("NOM_SD_%s", a), f4(d$sd_SC)); tok(sprintf("NOM_SEC_%s", a), f2(d$mean_secs)) }

# ---- E5. SENSITIVITY: population size under the equal budget (k = 3) --------------------
if (RUN_SENSITIVITY) {
  sens <- do.call(rbind, lapply(META, function(a) do.call(rbind, lapply(c(20, 50, 100), function(ps) {
    v <- sapply(seq_len(N_RUN), function(r) run_one(a, K_FINAL, SEED + r, pop = ps, target = GLOBAL[["3"]])$best)
    data.frame(algorithm = a, population = ps, mean_SC = mean(v), sd_SC = sd(v),
               success_rate = mean(v >= GLOBAL[["3"]] - 1e-9)) }))))
  wcsv(sens, "E5_sensitivity_population_k3.csv"); print(sens, row.names = FALSE)
  rngs <- sapply(META, function(a) diff(range(sens$mean_SC[sens$algorithm == a])))
  tok("SENS_SENTENCE", paste0("varying the population or swarm size over {20, 50, 100} at the same 3,000-evaluation budget changed the mean k = 3 Silhouette by at most ",
    paste(sprintf("%s (%s)", f4(rngs), names(rngs)), collapse = ", ")))
} else tok("SENS_SENTENCE", "sensitivity analysis not run")

# ---- E6. FINAL PARTITIONS (global optima) --------------------------------------------
lab3 <- labels_of(glob[["3"]]$med); lab2 <- labels_of(glob[["2"]]$med)
relabel3 <- function(lab) { cj <- lab[PU == "DKI JAKARTA"]; cb <- lab[PU == "BANTEN"]
  ifelse(lab == cj, 1L, ifelse(lab == cb, 3L, 2L)) }
L3 <- relabel3(lab3)
sw3 <- cluster::silhouette(lab3, DIST); sw2 <- cluster::silhouette(lab2, DIST)
memb <- data.frame(Province = provinsi, Cluster_k3 = L3, s_i_k3 = round(sw3[, 3], 4),
                   Cluster_k2 = lab2, s_i_k2 = round(sw2[, 3], 4))
wcsv(memb, "E6_membership_and_silhouette_k2_k3.csv")
prof <- aggregate(data_num, by = list(Cluster = L3), FUN = mean)
prof$Average_of_4 <- rowMeans(prof[, SHORT]); prof$n <- as.integer(table(L3))
wcsv(cbind(prof[1], round(prof[-1], 2)), "TABLE4_profile_k3.csv"); print(round(prof, 2))
tok("K3_SIZES", paste(as.integer(table(L3)), collapse = "\u2013"))
tok("K3_MEDOIDS", paste(provinsi[glob[["3"]]$med], collapse = ", "))
tok("K3_SIL_BIG", f4(mean(sw3[L3 == 2, 3])))
sz2 <- table(lab2); minor <- as.integer(names(sz2)[which.min(sz2)])
tok("K2_SIZES", paste(sort(as.integer(sz2)), collapse = "\u2013"))
tok("K2_MINOR", paste(provinsi[lab2 == minor], collapse = ", "))
tok("K3_NEG_SI", sum(sw3[, 3] < 0))

# ---- E7. DBI / CH ACROSS k (best-known partitions) ------------------------------------
XM <- as.matrix(X)
idx_k <- do.call(rbind, lapply(K_RANGE, function(k) {
  med <- if (as.character(k) %in% names(glob)) glob[[as.character(k)]]$med else best_part[[as.character(k)]]$med
  lab <- labels_of(med)
  cr <- clusterCrit::intCriteria(XM, as.integer(lab), c("Davies_Bouldin", "Calinski_Harabasz"))
  data.frame(k = k, SC = sil_med(med), DBI = cr$davies_bouldin, CH = cr$calinski_harabasz,
             sizes = paste(sort(tabulate(lab, k)), collapse = "-"),
             source = if (as.character(k) %in% names(glob)) "global optimum" else "best of equal-budget runs")
}))
wcsv(idx_k, "E7_indices_across_k.csv"); print(idx_k, row.names = FALSE)
tok("DBI_K3", f4(idx_k$DBI[idx_k$k == 3])); tok("CH_K3", f3(idx_k$CH[idx_k$k == 3]))
tok("DBI_ARGMIN_K", idx_k$k[which.min(idx_k$DBI)]); tok("CH_ARGMAX_K", idx_k$k[which.max(idx_k$CH)])
tok("SC_ARGMAX_K", idx_k$k[which.max(idx_k$SC)])

# ---- E8. LEAVE-ONE-PROVINCE-OUT STABILITY -------------------------------------------
if (RUN_LOO) {
  loo <- do.call(rbind, lapply(seq_len(n), function(j) {
    Xj <- scale(data_num[-j, ]); Dj <- as.matrix(dist(Xj)); nj <- n - 1
    res <- sapply(2:3, function(k) {
      cm <- utils::combn(nj, k)
      sc <- vapply(seq_len(ncol(cm)), function(c) fast_sil_lab(labels_of(cm[, c], Dj), Dj), numeric(1))
      labj <- labels_of(cm[, which.max(sc)], Dj)
      full <- if (k == 3) lab3[-j] else lab2[-j]
      c(ari(full, labj), max(sc))
    })
    data.frame(removed = provinsi[j], ARI_k2 = res[1, 1], ARI_k3 = res[1, 2], SC_k2 = res[2, 1], SC_k3 = res[2, 2])
  }))
  wcsv(loo, "E8_leave_one_out_stability.csv")
  ch3 <- loo$removed[loo$ARI_k3 < 1 - 1e-9]
  tok("LOO_K3_MEAN", f3(mean(loo$ARI_k3))); tok("LOO_K3_MIN", f3(min(loo$ARI_k3)))
  tok("LOO_K3_NSAME", sum(loo$ARI_k3 > 1 - 1e-9))
  tok("LOO_K3_CHANGED", if (length(ch3)) paste(ch3, collapse = ", ") else "none")
  tok("LOO_K2_MEAN", f3(mean(loo$ARI_k2))); tok("LOO_K2_NSAME", sum(loo$ARI_k2 > 1 - 1e-9))
}

# ---- E9. SILHOUETTE-vs-k FIGURE ------------------------------------------------------
sk <- rbind(data.frame(k = K_RANGE, method = "PAM baseline", SC = pam_tab$PAM_SC),
            data.frame(k = as.integer(names(GLOBAL)), method = "Global optimum (exhaustive)", SC = unname(GLOBAL)),
            do.call(rbind, lapply(ALGS_EB, function(a) data.frame(k = K_RANGE, method = paste0(a, " (mean, 3,000 evals)"),
                                                                  SC = sapply(K_RANGE, function(k) G(k, a, "mean_SC"))))))
sk$method <- factor(sk$method, levels = unique(sk$method))
psk <- ggplot(sk, aes(k, SC, colour = method, shape = method)) + geom_line() + geom_point(size = 2) +
  geom_vline(xintercept = K_FINAL, linetype = 3) + scale_x_continuous(breaks = K_RANGE) +
  scale_shape_manual(values = c(4, 8, 16, 17, 15, 18, 1)) +
  labs(x = "Number of clusters k", y = "Mean Silhouette Coefficient", colour = NULL, shape = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) + theme(legend.text = element_text(size = 8))
ggsave(file.path(DIR_OUT, "Fig_silhouette_vs_k.png"), psk, width = 7, height = 5, dpi = 300)

# ---- E10. MAP (strict name matching, Java inset) -------------------------------------
shp <- find_shp(PATH_SHP)
if (HAS_SF && !is.null(shp)) tryCatch({
  say("Map: using shapefile %s", shp)
  sfx <- sf::st_read(shp, quiet = TRUE)
  cand_cols <- names(sfx)[sapply(sf::st_drop_geometry(sfx), function(v) is.character(v) || is.factor(v))]
  norm <- function(x) { x <- toupper(trimws(gsub("\\s+", " ", x)))
    x <- gsub("-", " ", x); x <- gsub("^KEP\\.? ", "KEPULAUAN ", x)
    x <- gsub("^(KEPULAUAN )?BANGKA BELITUNG$", "KEPULAUAN BANGKA BELITUNG", x); x <- gsub("^YOGYAKARTA$", "DI YOGYAKARTA", x); x <- gsub("^(D\\.?I\\.?|DAERAH ISTIMEWA) YOGYAKARTA$", "DI YOGYAKARTA", x)
    x <- gsub("^DKI JAKARTA$|^DAERAH KHUSUS IBUKOTA JAKARTA$|^JAKARTA RAYA$", "DKI JAKARTA", x)
    gsub("^NANGGROE ACEH DARUSSALAM$", "ACEH", x) }
  score <- sapply(cand_cols, function(cn) sum(norm(sf::st_drop_geometry(sfx)[[cn]]) %in% norm(provinsi)))
  name_col <- cand_cols[which.max(score)]
  sfx$KEY <- norm(sf::st_drop_geometry(sfx)[[name_col]])
  key_bps <- norm(provinsi)
  parent <- c("PAPUA SELATAN" = "PAPUA", "PAPUA TENGAH" = "PAPUA", "PAPUA PEGUNUNGAN" = "PAPUA", "PAPUA BARAT DAYA" = "PAPUA BARAT")
  if (!all(names(parent) %in% sfx$KEY)) {                 # pre-2022 shapefile (34 provinces)
    kids <- key_bps %in% names(parent)
    chk_par <- tapply(L3, ifelse(kids, parent[key_bps], key_bps), function(v) length(unique(v)))
    if (any(chk_par > 1)) stop("Shapefile merges new Papua provinces that belong to different clusters.")
    key_bps <- ifelse(kids, parent[key_bps], key_bps)
  }
  cl <- data.frame(KEY = key_bps, Cluster = L3)[!duplicated(key_bps), ]
  unmatched_bps <- setdiff(cl$KEY, sfx$KEY); unmatched_shp <- setdiff(sfx$KEY, cl$KEY)
  wcsv(data.frame(side = c(rep("BPS not in shapefile", length(unmatched_bps)), rep("shapefile not in BPS", length(unmatched_shp))),
                  name = c(unmatched_bps, unmatched_shp)), "E10_map_name_matching_report.csv")
  if (length(unmatched_bps) || length(unmatched_shp))
    stop("Map: unmatched names (see E10_map_name_matching_report.csv). Fix recodes in norm(); no silent fallback.")
  sfx <- merge(sfx, cl, by = "KEY", all.x = TRUE)
  chk_map <- table(sf::st_drop_geometry(sfx)$Cluster)
  print(chk_map)
  sfx$Cluster <- factor(sfx$Cluster, levels = 1:3,
    labels = c(sprintf("Cluster 1: DKI Jakarta (mean of 4 components %.2f)", prof$Average_of_4[1]),
               sprintf("Cluster 2: %d provinces (%.2f)", prof$n[2], prof$Average_of_4[2]),
               sprintf("Cluster 3: Banten (%.2f)", prof$Average_of_4[3])))
  pal <- c("#d73027", "#4d9221", "#fdae61")
  base <- ggplot(sfx) + geom_sf(aes(fill = Cluster), colour = "grey30", linewidth = 0.1) +
    scale_fill_manual(values = pal, drop = FALSE, name = NULL)
  main <- base + coord_sf(expand = FALSE) + labs(x = "Longitude (\u00b0E)", y = "Latitude") +
    annotate("rect", xmin = 105, xmax = 107.8, ymin = -7.3, ymax = -5.7, fill = NA, colour = "black", linewidth = 0.4) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom", legend.direction = "vertical",
                                     legend.text = element_text(size = 16), legend.key.size = grid::unit(0.8, "cm"),
                                     panel.grid = element_line(colour = "grey90", linewidth = 0.2))
  inset <- base + coord_sf(xlim = c(105, 107.8), ylim = c(-7.3, -5.7), expand = FALSE) +
    annotate("label", x = c(106.83, 106.1), y = c(-5.95, -6.65), label = c("DKI Jakarta", "Banten"), size = 2.6) +
    theme_void() + theme(legend.position = "none", panel.border = element_rect(fill = NA, colour = "black"),
                         plot.background = element_rect(fill = "white", colour = NA))
  png(file.path(DIR_OUT, "Fig_cluster_map_k3.png"), width = 12, height = 6.2, units = "in", res = 300)
  print(main)
  print(inset, vp = grid::viewport(x = 0.83, y = 0.80, width = 0.28, height = 0.32))
  dev.off()
  tok("MAP_NPOLY", nrow(sfx))
}, error = function(e) { if (dev.cur() > 1) dev.off()
  message("MAP NOT PRODUCED (all other outputs are fine): ", conditionMessage(e)) }
) else message("Map skipped: no Indonesian province shapefile found. Set PATH_SHP to the 34/38-province .shp (not a Jawa Tengah file).")

# ---- SESSION, PLACEHOLDERS ----------------------------------------------------------
tok("R_VERSION", paste(R.version$major, R.version$minor, sep = ".")); tok("R_YEAR", R.version$year)
for (pk in c("GA", "pso", "cluster", "clusterCrit")) tok(paste0("V_", toupper(pk)), as.character(packageVersion(pk)))
tok("PLATFORM", paste(Sys.info()[["sysname"]], R.version$platform))
tok("N_RUN", N_RUN); tok("BUDGET", fint(BUDGET))
tok("TOTAL_MIN", sprintf("%.0f", (proc.time()[["elapsed"]] - T0_ALL) / 60))
ph <- data.frame(key = ls(TOK), value = unlist(mget(ls(TOK), envir = TOK)))
write.csv(ph, file.path(DIR_OUT, "placeholders.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(capture.output(sessionInfo()), file.path(DIR_OUT, "R_sessionInfo.txt"))
say("\nDONE in %.1f min. Outputs: %s", (proc.time()[["elapsed"]] - T0_ALL) / 60, normalizePath(DIR_OUT))
say("Upload placeholders.csv + the Fig_*.png files to fill the revised manuscript automatically.")
