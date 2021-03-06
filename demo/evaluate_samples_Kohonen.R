# satellite image time series package (SITS)
# example of clustering using self-organizin maps
library(sits)

# load the sitsdata library
if (!requireNamespace("sitsdata", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools")
    }
    devtools::install_github("e-sensing/sitsdata")
}
library(sitsdata)
data("samples_cerrado_mod13q1")

# Clustering time series samples using self-organizing maps
som_map <-
    sits_som_map(
        samples_cerrado_mod13q1,
        grid_xdim = 12,
        grid_ydim = 12,
        alpha = 1,
        distance = "euclidean"
    )

plot(som_map)

# Remove samples that have  bad quality
new_samples <- sits_som_clean_samples(som_map)

cluster_purity <- sits_som_evaluate_cluster(som_map)
plot(cluster_purity)
