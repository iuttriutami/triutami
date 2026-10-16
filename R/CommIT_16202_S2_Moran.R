# =============================================================================
# CommIT 16202 - SUPPLEMENTARY S2: Global Moran's I of the EQI 2023 components
# (Reviewer 3: the map alone does not establish spatial autocorrelation)
# Weights: k-nearest neighbours on province locations (k = 4 main; 3, 5, 6 sensitivity),
# row-standardised. p-values: analytic (randomisation) and permutation (9,999).
# Province location = mean of the centroids of its regencies/cities under the 2024
# administrative codes (PROV_2024), so all 38 provinces, including the four Papua
# provinces created in 2022, are represented. No shapefile is needed.
# =============================================================================
PATH_CSV <- "/kaggle/input/datasets/triutamiiut/2023-eqi/Komponen Penyusun Indeks Kualitas Lingkungan Hidup Menurut Provinsi 2023.csv"
PATH_XY  <- "/kaggle/input/datasets/triutamiiut/xy-centroid/xy centroid.csv"
SEED <- 42L; KNN_K <- 3:6; NSIM <- 9999
DIR_OUT <- if (dir.exists("/kaggle/working")) "/kaggle/working/commIT_supplementary" else "commIT_supplementary"
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

for (p in c("spdep")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
find_in <- function(path, pattern) {
  if (file.exists(path)) return(path)
  cand <- list.files("/kaggle/input", pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (!length(cand)) stop("File not found (", pattern, "): add the dataset via 'Add Input' and rerun.")
  message("Using: ", cand[1]); cand[1]
}
PATH_CSV <- find_in(PATH_CSV, "Kualitas Lingkungan.*\\.csv$")
# Province coordinates: either the regency/city centroid file (PATH_XY, columns PROV_2024,
# X_Centroid, Y_Centroid), or, if it is not available, the province coordinates derived from it
# and distributed with the repository (results/supplementary/supp_S2_coordinates.csv).
PATH_PROV_XY <- "results/supplementary/supp_S2_coordinates.csv"
USE_PROV_XY  <- !file.exists(PATH_XY) && !length(list.files("/kaggle/input", pattern = "centroid.*\\.csv$", recursive = TRUE)) &&
                file.exists(PATH_PROV_XY)
if (!USE_PROV_XY) PATH_XY <- find_in(PATH_XY, "centroid.*\\.csv$")

# ---- data (same four components as the main analysis) ----
df <- read.csv(PATH_CSV, header = TRUE, check.names = FALSE, na.strings = c("", "NA", "-"))
names(df)[1] <- "Provinsi"
df <- df[!toupper(trimws(df$Provinsi)) %in% c("INDONESIA", ""), , drop = FALSE]
nm <- names(df)
pick <- function(ok) { j <- which(ok); if (length(j) != 1) stop("Column not identified: ", paste(nm[j], collapse = "; ")); j }
j <- c(AQI = pick(grepl("Kualitas Udara", nm)), WQI = pick(grepl("Kualitas Air", nm) & !grepl("Laut", nm)),
       LCQI = pick(grepl("Tutupan Lahan", nm)), LQI = pick(grepl("Kualitas Lahan", nm)))
X <- as.data.frame(lapply(df[, j], function(v) as.numeric(gsub(",", ".", as.character(v))))); names(X) <- names(j)
prov <- toupper(trimws(gsub("\\s+", " ", as.character(df$Provinsi))))
stopifnot(length(prov) == 38L, !anyNA(X))
vars <- cbind(X, Average_of_4 = rowMeans(X))

# ---- locations: mean of regency/city centroids per province (2024 codes, 38 provinces) ----
if (USE_PROV_XY) {
  message("Using province coordinates from ", PATH_PROV_XY)
  pxy <- read.csv(PATH_PROV_XY, check.names = FALSE)
  norm <- function(x) toupper(trimws(gsub("\\s+", " ", x)))
  key <- norm(prov)
  xy  <- as.matrix(pxy[match(key, norm(pxy$Province)), c("lon", "lat")]); rownames(xy) <- key
} else {
  cen <- read.csv(PATH_XY, check.names = FALSE)
  norm <- function(x) { x <- toupper(trimws(gsub("\\s+", " ", x))); x <- gsub("-", " ", x)
    x <- gsub("^KEP\\.? ", "KEPULAUAN ", x); x <- gsub("^(KEPULAUAN )?BANGKA BELITUNG$", "KEPULAUAN BANGKA BELITUNG", x)
    x <- gsub("^YOGYAKARTA$|^(D\\.?I\\.?|DAERAH ISTIMEWA) YOGYAKARTA$", "DI YOGYAKARTA", x)
    x <- gsub("^DAERAH KHUSUS IBUKOTA JAKARTA$|^JAKARTA RAYA$|^D\\.?K\\.?I\\.? JAKARTA$", "DKI JAKARTA", x)
    gsub("^NANGGROE ACEH DARUSSALAM$", "ACEH", x) }
  cen$KEY <- norm(cen$PROV_2024)
  pc <- aggregate(cbind(lon = X_Centroid, lat = Y_Centroid) ~ KEY, data = cen, FUN = mean)
  key <- norm(prov)
  miss <- setdiff(key, pc$KEY)
  if (length(miss)) { print(sort(unique(pc$KEY))); stop("Not in centroid file: ", paste(miss, collapse = ", ")) }
  xy <- as.matrix(pc[match(key, pc$KEY), c("lon", "lat")]); rownames(xy) <- key
}
stopifnot(!anyNA(xy), nrow(xy) == 38L)
if (!USE_PROV_XY) write.csv(data.frame(Province = prov, lon = xy[, 1], lat = xy[, 2],
                     n_regencies = as.integer(table(cen$KEY)[key])),
          file.path(DIR_OUT, "supp_S2_coordinates.csv"), row.names = FALSE)

# ---- Moran's I ----
mor <- do.call(rbind, lapply(KNN_K, function(kk) {
  lw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(xy, k = kk, longlat = TRUE)), style = "W")
  do.call(rbind, lapply(names(vars), function(v) {
    mt <- spdep::moran.test(vars[[v]], lw, randomisation = TRUE)
    set.seed(SEED); mc <- spdep::moran.mc(vars[[v]], lw, nsim = NSIM)
    data.frame(knn = kk, variable = v, I = unname(mt$estimate[1]), E_I = unname(mt$estimate[2]),
               z = unname(mt$statistic), p_randomisation = mt$p.value, p_permutation = mc$p.value)
  }))
}))
write.csv(transform(mor, I = round(I, 4), E_I = round(E_I, 4), z = round(z, 3)),
          file.path(DIR_OUT, "supp_S2_moran_2023.csv"), row.names = FALSE)
cat("\n=== Moran's I, k = 4 nearest neighbours (main) ===\n")
print(transform(mor[mor$knn == 4, -1], I = round(I, 3), z = round(z, 2)), row.names = FALSE)
cat("\n=== Significant at 0.05 (permutation p) for each k ===\n")
print(tapply(mor$p_permutation < 0.05, list(mor$variable, mor$knn), identity))
writeLines(capture.output(sessionInfo()), file.path(DIR_OUT, "supp_S2_sessionInfo.txt"))
cat("\nDONE ->", normalizePath(DIR_OUT), "\n")
