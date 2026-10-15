# Metaheuristic-Optimized K-Medoids Clustering of Indonesian Provinces (EQI 2023)

Code, data and results for the manuscript *“Metaheuristic-Optimized K-Medoids Clustering of Indonesian Provinces Based on Environmental Quality Index Components”* (CommIT Journal, article 16202).

![R Shiny decision-support application](docs/shiny_app_screenshot.png)

## Contents

| Folder | What it contains |
|---|---|
| `R/CommIT_16202_Revision_Experiments.R` | Main pipeline: preprocessing (Mahalanobis D², VIF), PAM baseline, exhaustive global optimum (k = 2–6), equal-budget benchmark of GA, PSO, CSA, GWO and random search (3,000 fitness evaluations, 30 runs, k = 2–10), nominal-budget runs, sensitivity analysis, leave-one-province-out stability (ARI), DBI/CH, figures and map |
| `R/CommIT_16202_S1_EQI2022.R` | Supplementary S1: replication of the full protocol on EQI 2022 (34 provinces) |
| `R/CommIT_16202_S2_Moran.R` | Supplementary S2: global Moran’s I of the 2023 components (k-nearest-neighbour weights, 9,999 permutations) |
| `data/` | BPS tables *Komponen Penyusun Indeks Kualitas Lingkungan Hidup Menurut Provinsi* (2022 and 2023) |
| `results/main/` | CSV outputs of the main pipeline (E0–E8, Tables 2 and 4) |
| `results/supplementary/` | CSV outputs of S1 (EQI 2022) and S2 (Moran’s I) |
| `figures/` | Figures used in the manuscript (PNG) |
| `shiny/app.R` | Interactive R Shiny decision-support application (English interface) |

## How to reproduce

1. Install R ≥ 4.4 and the packages `cluster`, `GA`, `pso`, `clusterCrit`, `ggplot2`, `sf`, `spdep`.
2. Edit the paths at the top of each script (`PATH_CSV`, `PATH_SHP`, `PATH_XY`) to point to `data/`, or run the scripts on Kaggle with the data under `/kaggle/input`.
3. Run `R/CommIT_16202_Revision_Experiments.R` (about 40 minutes on a Kaggle CPU notebook), then the S1 and S2 scripts.

Every run of every optimizer uses seed 42 + r (r = 1, …, 30). The fitness is the mean Silhouette with nearest-medoid assignment; its vectorized implementation is checked against `cluster::silhouette` before any experiment runs. Silhouette values are rounded to 10 decimal places before the Friedman and Wilcoxon tests, so that equivalent medoid sets are treated as ties.

## Run the Shiny app locally

```r
install.packages(c("shiny", "shinydashboard", "DT", "ggplot2", "cluster", "GA", "pso", "clusterCrit", "dplyr", "tidyr"))
shiny::runApp("shiny")
```

Online version: https://triutami.shinyapps.io/myshinyapp/

## Data source

Badan Pusat Statistik (Statistics Indonesia), *Komponen Penyusun Indeks Kualitas Lingkungan Hidup Menurut Provinsi*, 2022 and 2023, https://www.bps.go.id/id/statistics-table/2/MjUzMiMy/komponen-penyusun-indeks-kualitas-lingkungan-hidup-menurut-provinsi.html

## Citation

Please cite the article once published (CommIT Journal, article ID 16202).
