# =============================================================================
# Check that every optimizer in shiny/app.R returns the medoid set that produced
# the reported Silhouette (no stochastic re-decoding).
# Usage (from the repository root):  Rscript shiny/check_shiny_decoding.R
# =============================================================================
suppressPackageStartupMessages({
  for (p in c("shiny", "shinydashboard", "DT", "ggplot2", "dplyr", "tidyr", "cluster", "clusterCrit", "GA", "pso"))
    library(p, character.only = TRUE)
})
app_file <- if (file.exists("shiny/app.R")) "shiny/app.R" else "app.R"
exprs <- parse(app_file)
keep  <- !vapply(exprs, function(e) is.call(e) && identical(e[[1]], as.name("shinyApp")), logical(1))
for (e in exprs[keep]) eval(e, globalenv())          # load helpers and optimizers, do not start the app

csv <- list.files(c("data", "../data", "."), pattern = "Provinsi.*2023\\.csv$", full.names = TRUE)[1]
if (is.na(csv)) { if (interactive()) csv <- file.choose() else stop("2023 CSV not found in data/") }
skip <- grep("Kualitas Udara", readLines(csv, warn = FALSE))[1] - 1
df  <- read.csv(csv, skip = skip, check.names = FALSE, na.strings = c("", "NA", "-"))
names(df)[1] <- "Provinsi"
df  <- df[!toupper(trimws(df$Provinsi)) %in% c("INDONESIA", "") & !is.na(df$Provinsi), ]
nm  <- names(df)
j   <- c(which(grepl("Kualitas Udara", nm)), which(grepl("Kualitas Air", nm) & !grepl("Laut", nm)),
         which(grepl("Tutupan Lahan", nm)), which(grepl("Kualitas Lahan", nm)))
X   <- as.data.frame(lapply(df[, j], function(v) as.numeric(gsub(",", ".", as.character(v)))))
stopifnot(ncol(X) == 4, nrow(X) == 38, !anyNA(X))
S    <- scale(X); DIST <- dist(S); DMAT <- as.matrix(DIST); n <- nrow(S)

out <- list()
for (k in c(3, 5)) for (seed in 43:47) {
  for (a in c("GA", "PSO", "CSA", "GWO")) {
    r <- switch(a,
      GA  = run_GA (k, n, DMAT, DIST, seed = seed, popSize = 30, maxiter = 50),
      PSO = run_PSO(k, n, DMAT, DIST, seed = seed, s = 30, maxit = 50),
      CSA = run_CSA(k, n, DMAT, DIST, seed = seed, flock = 30, num_iter = 50),
      GWO = run_GWO(k, n, DMAT, DIST, seed = seed, num_wolves = 30, num_iter = 50))
    recomputed <- mean(silhouette(labels_of_medoids(DMAT, r$med), DIST)[, 3])
    out[[length(out) + 1]] <- data.frame(k = k, seed = seed, algorithm = a, reported = r$sil,
                                         recomputed = recomputed, abs_diff = abs(r$sil - recomputed))
  }
}
res <- do.call(rbind, out)
print(res, row.names = FALSE, digits = 6)
max_diff <- max(res$abs_diff)
status   <- if (max_diff < 1e-9) "PASS" else "FAIL"
cat(sprintf("\nMax |reported - recomputed| = %.2e  ->  %s\n", max_diff, status))

# ---- log (written next to this script) ----
log_dir <- dirname(normalizePath(app_file))
write.csv(res, file.path(log_dir, "check_shiny_decoding_log.csv"), row.names = FALSE)
writeLines(c(sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             sprintf("Data: %s", basename(csv)),
             sprintf("Runs: %d (k = 3 and 5; seeds 43-47; GA, PSO, CSA, GWO)", nrow(res)),
             sprintf("Max |reported - recomputed| = %.3e", max_diff),
             sprintf("Result: %s", status), "", capture.output(sessionInfo())),
           file.path(log_dir, "check_shiny_decoding_log.txt"))
stopifnot("Shiny decoding check failed: returned medoids do not reproduce the reported Silhouette" = max_diff < 1e-9)
