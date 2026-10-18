# Subgingival Ecological Landscape

This repository contains MATLAB and R code for investigating the ecological landscape of the subgingival microbiome and its associations with host characteristics.

The computational workflows include:

- Unsupervised feature selection using iDetect.
- Nonlinear manifold reconstruction using DDRTree principal-tree learning.
- Spectral clustering using distances on the learned tree.
- Ridge regularized canonical correlation analysis (rCCA) relating host metadata to bacterial taxa or microbiome cluster membership.

The landscape analysis was developed to investigate discrete microbiome states within continuous ecological transitions from eubiosis to dysbiosis in a population-based cohort.

## Repository contents

| Path | Description |
|---|---|
| `demo_run.m` | MATLAB entry point. |
| `run_landscape.m` | Main landscape analysis function. |
| `config_example.m` | Input paths, output paths, model parameters, and figure settings. |
| `data/MIP_data_ready.mat` | Prepared MATLAB input included in this repository snapshot. |
| `dependencies/` | iDetect, DDRTree, and landscape plotting implementations. |
| `results/` | MATLAB model outputs, cluster assignments, feature weights, and figures. |
| `rCCA/rCCA_Species.R` | R Markdown-formatted analysis of host metadata and bacterial taxa. |
| `rCCA/rCCA_Clusters.R` | R Markdown-formatted analysis of cluster membership and host metadata. |
| `rCCA/Data/` | Small previews of metadata, bacterial values, and cluster coding. |
| `rCCA/Results/` | Stored rCCA loading and variate tables. |
| `rCCA/sessionInfo.txt` | R version, platform, and package versions for the recorded analysis session. |

## MATLAB landscape analysis

### Input

The analysis expects a MATLAB MAT-file containing a structure named `x_ready`:

| Field | Description |
|---|---|
| `x_ready.clr` | Numeric taxa-by-sample matrix of centered log-ratio (CLR)-transformed abundances. |
| `x_ready.sample_id` | Sample identifiers corresponding to matrix columns. |
| `x_ready.tax` | Taxon names corresponding to matrix rows; recommended but optional. |

The default input is `data/MIP_data_ready.mat`. Change `config.input_mat_file` in `config_example.m` to use another file.

### Dependencies

- MATLAB with the Statistics and Machine Learning Toolbox.
- `dependencies/MIP_iDetect.m`.
- `dependencies/MIP_DDRTree.m`.
- `dependencies/plotDDRtree.m`.

`MIP_DDRTree.m` also calls `pca_projection` and `sqdist`. These helper functions are not included in this repository snapshot and must be available on the MATLAB path.

### Running the analysis

Set the MATLAB current folder to the repository root, review `config_example.m`, and run:

```matlab
config = config_example();
results = run_landscape(config);
```

The supplied configuration uses a three-dimensional embedding and eight clusters. Input paths, output locations, and model parameters are defined in `config_example.m`.

Unsupervised cluster labels are arbitrary. Verify `config.clustering.cluster_order` after regenerating the model before interpreting the reordered labels.

### Outputs

- `01_idetect_feature_selection.mat`
- `02_ddrtree_model.mat`
- `03_tree_geodesic_distance.mat`
- `MIP_subgingival_landscape_results.mat`
- `sample_cluster_assignments.tsv`
- `idetect_feature_weights.tsv`
- `idetect_feature_weights.pdf`
- Landscape figures in PDF and MATLAB figure formats.

Outputs are written to the configured directory. Running with the default settings can replace files with matching names in `results/`.

## R regularized canonical correlation analysis

### Analyses

Both analyses use `mixOmics::rcc()` with ridge regularization.

| Analysis | X block | Y block |
|---|---|---|
| Species | Cleaned and imputed host metadata | CLR-transformed bacterial abundances |
| Clusters | Cluster membership coding | Cleaned and imputed host metadata |

The files include canonical correlation plots, variable plots, canonical variates, correlation-based loadings, clustered image maps, and a five-fold regularization grid search. Some sections remain placeholders and require preparation before the documents can run from beginning to end.

### R environment

The recorded session used **R 4.6.1 on Windows 11 x64**. Full session information, including attached packages and loaded namespaces, is provided in [`rCCA/sessionInfo.txt`](rCCA/sessionInfo.txt).

Selected package versions from that session are listed below:

| Package | Version |
|---|---|
| mixOmics | 6.36.0 |
| vegan | 2.7-6 |
| ggplot2 | 4.0.3 |
| dplyr | 1.2.1 |
| igraph | 2.3.4 |
| plotly | 4.12.1 |
| readxl | 1.5.0.1 |
| openxlsx | 4.2.9 |
| writexl | 2.0.1 |

The session record documents the environment; it does not automatically install or restore these versions.

### Input preparation

The species analysis requires numeric sample-by-variable tables of host metadata and CLR-transformed bacterial abundances.

The cluster analysis requires host metadata and numeric cluster membership coding, with one indicator column per cluster (one-hot encoded).

The small tables in `rCCA/Data/` illustrate the format...


### Running the R analyses

1. Save copies of `rCCA_Species.R` and `rCCA_Clusters.R` with `.Rmd` extensions. Their contents use R Markdown headers and code chunks.
2. Open the `.Rmd` files in RStudio and install the packages listed in their setup chunks. The species file includes installation of `mixOmics` through `BiocManager`.
3. Replace each placeholder with an input path. Ensure each reader matches the actual file format. 


### Stored results

| File in `rCCA/Results/` | Description |
|---|---|
| `loadings_X.xlsx` | Species analysis: X-block loading table. |
| `loadings_Y.xlsx` | Species analysis: Y-block loading table. |
| `variates_XY.xlsx` | Species analysis: stored variate table. |
| `clusters_loadings_X.xlsx` | Cluster analysis: X-block loading table. |
| `clusters_loadings_Y.xlsx` | Cluster analysis: Y-block loading table. |
| `clusters_variates_XY.xlsx` | Cluster analysis: stored variate table. |


## Data availability

This repository snapshot includes a prepared MATLAB input file, R input previews, and stored analysis outputs.

## License

See [`LICENSE`](LICENSE) for the repository's MIT license. Third-party code and dependencies remain subject to their applicable licenses.
