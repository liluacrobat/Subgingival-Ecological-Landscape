# Subgingival Ecological Landscape
This repository contains the MATLAB code used to reconstruct the ecological landscape of the subgingival microbiome using nonlinear manifold learning and graph-based clustering.
The computational framework integrates:  
	•	Unsupervised feature selection  
	•	Nonlinear manifold reconstruction (DDRTree principal-tree learning)  
	•	Spectral clustering  
	
The pipeline was developed to identify discrete microbiome states embedded within continuous ecological transitions from eubiosis to dysbiosis in a large population-based cohort.

## Repository structure

```text
.
├── run_landscape.m       Main analysis function
├── config_example.m      User-editable configuration
├── data/                 Input data directory; data are not distributed
├── dependencies/         External or study-specific MATLAB functions
└── results/              Generated outputs
```

## Input

The analysis expects a MATLAB MAT-file containing a structure named `x_ready` with the following fields:

- `x_ready.clr`: taxa-by-sample centered log-ratio transformed abundance matrix
- `x_ready.sample_id`: sample identifiers
- `x_ready.tax`: taxon names, recommended but optional

Participant-level data are not included in this repository. Place the prepared input file in `data/` or update `config.input_mat_file` to its local location.

## Dependencies

The following functions or toolboxes are required:

- MATLAB Statistics and Machine Learning Toolbox
- `MIP_iDetect.m`
- `MIP_DDRTree.m`
- `plotDDRtree.m`

Place the required study-specific implementations and their dependencies in the `dependencies/` directory. The repository should include only code that can be redistributed under its applicable license.

## Running the analysis

From MATLAB:

```matlab
config = config_example();
results = run_landscape(config);
```

All paths and model parameters are defined in `config_example.m`. Outputs are written to the configured `results/` directory.

## Principal outputs

- `01_idetect_feature_selection.mat`
- `02_ddrtree_model.mat`
- `03_tree_geodesic_distance.mat`
- `MIP_subgingival_landscape_results.mat`
- `sample_cluster_assignments.tsv`
- `idetect_feature_weights.tsv`
- `idetect_feature_weights.pdf`
- `subgingival_landscape.pdf`
