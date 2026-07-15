function config = config_example()
%CONFIG_EXAMPLE Configuration for the MIP landscape analysis.
%
% Edit the paths below before running:
%   config = config_example();
%   results = run_landscape(config);

repo_dir = fileparts(mfilename('fullpath'));

config.input_mat_file = fullfile(repo_dir, 'data', 'MIP_data_ready.mat');
config.input_variable = 'x_ready';
config.dependency_dir = fullfile(repo_dir, 'dependencies');
config.output_dir = fullfile(repo_dir, 'results');

config.minimum_selected_features = 2;

% iDetect parameters used in the final MIP landscape analysis.
config.idetect.max_iterations = 100;
config.idetect.distance = 'block';
config.idetect.soft = 0;
config.idetect.alpha = 1;
config.idetect.kernel = 'parabolic';
config.idetect.sigma = 10;
config.idetect.lambda = 15;
config.idetect.weight_threshold = 0.005;

% DDRTree parameters used in the final MIP landscape analysis.
config.ddrtree.sigma = 50;
config.ddrtree.lambda = 300;
config.ddrtree.gamma = 1;
config.ddrtree.max_iterations = 100;
config.ddrtree.relative_tolerance = 1e-9;
config.ddrtree.dimension = 3;

% Spectral clustering parameters.
config.clustering.similarity_sigma = 1;
config.clustering.number_of_clusters = 8;

% Original spectral-cluster labels reordered into the study's final
% landscape sequence. Verify this mapping after regenerating the model,
% because unsupervised cluster labels are arbitrary.
config.clustering.cluster_order = [4 5 2 3 1 8 7 6];

% Figure settings.
config.figures.create = true;
config.figures.plot_box_aspect_ratio = [1 2 1];
config.figures.azimuth = -83;
config.figures.elevation = 11;
end
