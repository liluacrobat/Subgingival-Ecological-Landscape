function results = run_landscape(config)
%RUN_LANDSCAPE Build the MIP subgingival microbiome landscape.
%
%   RESULTS = RUN_LANDSCAPE(CONFIG) performs the analysis used to construct
%   the subgingival microbiome landscape from a centered log-ratio (CLR)
%   transformed microbial abundance matrix. The workflow consists of:
%
%     1. iDetect-based feature selection
%     2. DDRTree principal-tree learning
%     3. Geodesic distance calculation along the learned tree
%     4. Spectral clustering of samples using tree-based distances
%
%   CONFIG is a structure returned by config_example.m. At minimum, the
%   input MAT-file must contain a structure named "x_ready" with fields:
%
%       x_ready.clr        Taxa-by-sample CLR-transformed abundance matrix
%       x_ready.tax        Taxon names (optional but recommended)
%       x_ready.sample_id  Sample identifiers
%
%   The function saves intermediate and final results to CONFIG.output_dir.
%
%   This code reproduces the landscape-building component of the MIP study.
%   It depends on the study implementation of iDetect and DDRTree, together
%   with MATLAB Statistics and Machine Learning Toolbox.
%
%   Author: Lu Li
%   Repository version: 1.0

arguments
    config (1,1) struct
end

config = validate_config(config);
addpath(genpath(config.dependency_dir));
if ~exist(config.output_dir, 'dir')
    mkdir(config.output_dir);
end

%% Load input data
input_data = load(config.input_mat_file, config.input_variable);
if ~isfield(input_data, config.input_variable)
    error('Input file does not contain variable "%s".', config.input_variable);
end
x_ready = input_data.(config.input_variable);
validate_input_data(x_ready);

clr_abundance = double(x_ready.clr);
if any(~isfinite(clr_abundance), 'all')
    error('The CLR abundance matrix contains NaN or Inf values.');
end

%% Stage 1: iDetect feature selection
fprintf('Stage 1/4: iDetect feature selection...\n');
idetect_parameters = struct( ...
    'it', config.idetect.max_iterations, ...
    'distance', config.idetect.distance, ...
    'soft', config.idetect.soft, ...
    'alpha', config.idetect.alpha, ...
    'kernel', config.idetect.kernel, ...
    'sigma', config.idetect.sigma, ...
    'lambda', config.idetect.lambda);
rng default
[objective, feature_weight, idetect_auxiliary] = ...
    MIP_iDetect(clr_abundance, idetect_parameters);

feature_weight = feature_weight(:);
if max(feature_weight) <= 0
    error('iDetect returned non-positive feature weights.');
end
normalized_weight = feature_weight ./ max(feature_weight);
selected_taxa = normalized_weight > config.idetect.weight_threshold;

if nnz(selected_taxa) < config.minimum_selected_features
    error(['Only %d features passed the iDetect threshold. Reduce the ', ...
        'threshold or review the input data.'], nnz(selected_taxa));
end

selected_clr = clr_abundance(selected_taxa, :);

idetect_results = struct;
idetect_results.objective = objective;
idetect_results.raw_weight = feature_weight;
idetect_results.normalized_weight = normalized_weight;
idetect_results.selected_taxa = selected_taxa;
idetect_results.auxiliary = idetect_auxiliary;
idetect_results.parameters = idetect_parameters;

save(fullfile(config.output_dir, '01_idetect_feature_selection.mat'), ...
    'idetect_results', '-v7.3');


%% Stage 2: DDRTree principal-tree learning
fprintf('Stage 2/4: DDRTree principal-tree learning...\n');
ddr_parameters = struct( ...
    'sigma', config.ddrtree.sigma, ...
    'lambda', config.ddrtree.lambda, ...
    'gamma', config.ddrtree.gamma, ...
    'maxIter', config.ddrtree.max_iterations, ...
    'eps', config.ddrtree.relative_tolerance, ...
    'dim', config.ddrtree.dimension);

[W, projection, edges, principal_tree, history] = ...
    MIP_DDRTree(selected_clr, ddr_parameters);

ddrtree_results = struct;
ddrtree_results.W = W;
ddrtree_results.projection = projection;
ddrtree_results.edges = edges;
ddrtree_results.principal_tree = principal_tree;
ddrtree_results.history = history;
ddrtree_results.parameters = ddr_parameters;

save(fullfile(config.output_dir, '02_ddrtree_model.mat'), ...
    'ddrtree_results', '-v7.3');


projection = ddrtree_results.projection;
edges = ddrtree_results.edges;
principal_tree = ddrtree_results.principal_tree;


%% Stage 3: Geodesic distances on the learned tree
fprintf('Stage 3/4: Computing tree geodesic distances...\n');
geodesic_distance = compute_tree_geodesic_distances( ...
    projection, principal_tree, edges);

save(fullfile(config.output_dir, '03_tree_geodesic_distance.mat'), ...
    'geodesic_distance', '-v7.3');

%% Stage 4: Spectral clustering
fprintf('Stage 4/4: Spectral clustering...\n');
similarity = exp(-(geodesic_distance ./ config.clustering.similarity_sigma).^2);


cluster_labels_raw = spectralcluster( ...
    similarity, config.clustering.number_of_clusters, ...
    'Distance', 'precomputed');

cluster_labels = reorder_cluster_labels( ...
    cluster_labels_raw, config.clustering.cluster_order);

%% Assemble and save final results
results = struct;
results.sample_id = x_ready.sample_id;
results.taxonomy = get_optional_field(x_ready, 'tax', strings(size(clr_abundance,1),1));
results.selected_taxa = selected_taxa;
results.selected_taxonomy = results.taxonomy(selected_taxa);
results.feature_weight = normalized_weight;
results.projection = projection;
results.principal_tree = principal_tree;
results.edges = edges;
results.geodesic_distance = geodesic_distance;
results.similarity = similarity;
results.cluster_labels_raw = cluster_labels_raw;
results.cluster_labels = cluster_labels;
results.config = config;

save(fullfile(config.output_dir, 'MIP_subgingival_landscape_results.mat'), ...
    'results', '-v7.3');

write_output_tables(results, config.output_dir);
create_summary_figures(results, config);

fprintf('Analysis complete. Results saved to: %s\n', config.output_dir);
end

function config = validate_config(config)
required = {'input_mat_file','input_variable','dependency_dir','output_dir', ...
    'minimum_selected_features','idetect','ddrtree','clustering'};
for i = 1:numel(required)
    if ~isfield(config, required{i})
        error('Missing required configuration field: config.%s', required{i});
    end
end

if ~isfile(config.input_mat_file)
    error('Input MAT-file not found: %s', config.input_mat_file);
end
if ~isfolder(config.dependency_dir)
    error('Dependency directory not found: %s', config.dependency_dir);
end
if config.clustering.number_of_clusters < 2
    error('The number of clusters must be at least 2.');
end
if numel(config.clustering.cluster_order) ~= config.clustering.number_of_clusters
    error('cluster_order must contain one entry per cluster.');
end
if ~isequal(sort(config.clustering.cluster_order(:))', ...
        1:config.clustering.number_of_clusters)
    error('cluster_order must be a permutation of 1:number_of_clusters.');
end
end

function validate_input_data(x_ready)
required = {'clr','sample_id'};
for i = 1:numel(required)
    if ~isfield(x_ready, required{i})
        error('Input structure x_ready is missing field "%s".', required{i});
    end
end
if size(x_ready.clr, 2) ~= numel(x_ready.sample_id)
    error('Number of CLR matrix columns must equal number of sample IDs.');
end
end

function geodesic_distance = compute_tree_geodesic_distances( ...
        projection, principal_tree, edges)
% Assign each sample to its nearest principal-tree node, then calculate
% pairwise shortest-path distances between assigned nodes.

edge_length_matrix = squareform(pdist(principal_tree')).^2 .* edges;
graph_object = graph(edge_length_matrix, 'upper');
node_distance = distances(graph_object);

nearest_node = knnsearch(principal_tree', projection');
geodesic_distance = node_distance(nearest_node, nearest_node);
geodesic_distance = (geodesic_distance+geodesic_distance')/2;
if any(~isfinite(geodesic_distance), 'all')
    error('The learned principal tree is disconnected.');
end
end

function reordered = reorder_cluster_labels(labels, order)
% Map the original cluster numbering to the prespecified study order.
reordered = zeros(size(labels));
for new_label = 1:numel(order)
    reordered(labels == order(new_label)) = new_label;
end
end

function value = get_optional_field(input_struct, field_name, default_value)
if isfield(input_struct, field_name)
    value = input_struct.(field_name);
else
    value = default_value;
end
end

function write_output_tables(results, output_dir)
sample_table = table(string(results.sample_id(:)), results.cluster_labels(:), ...
    'VariableNames', {'SampleID','Cluster'});
writetable(sample_table, fullfile(output_dir, 'sample_cluster_assignments.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

feature_table = table(string(results.taxonomy(:)), results.feature_weight(:), ...
    results.selected_taxa(:), ...
    'VariableNames', {'Taxon','NormalizedWeight','Selected'});
feature_table = sortrows(feature_table, 'NormalizedWeight', 'descend');
writetable(feature_table, fullfile(output_dir, 'idetect_feature_weights.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
end

function create_summary_figures(results, config)
if ~config.figures.create
    return
end

figure('Color','w');
positive_weights = results.feature_weight(results.feature_weight > 0);
plot(sort(positive_weights, 'descend'), '-k', 'LineWidth', 1.2);
hold on
x_limits = xlim;
plot(x_limits, repmat(config.idetect.weight_threshold, 1, 2), '--k');
set(gca, 'YScale', 'log', 'FontSize', 12);
xlabel('Features ranked by iDetect weight');
ylabel('Normalized feature weight');
title('iDetect feature selection');
exportgraphics(gcf, fullfile(config.output_dir, 'idetect_feature_weights.pdf'), ...
    'ContentType', 'vector');
close(gcf);

facecolorGroup = [0	0.776470588	0
0.003921569	0.098039216	0.576470588
0.015686275	0.2	1
0.741176471	0.215686275	0.91372549
1	0.580392157	0
1	0.149019608	0
0.650980392	0.466666667	0.156862745
0.5	0.5	0.5];
for i=1:max(results.cluster_labels)
    header{i} = ['Cluster ' num2str(i)];
end
plotDDRtree(results.projection, results.principal_tree, results.edges, results.cluster_labels, header,facecolorGroup);
title('MIP subgingival microbiome landscape');
pbaspect([1 2 1]);
axis([10 70 -50.0000 40.0000 -30.0000 30.0000]);
view(-83,11)
set(gca,'FontSize',18);
savefig(gcf, fullfile(config.output_dir, 'subgingival_landscape.fig'));
keyboard
exportgraphics(gcf, fullfile(config.output_dir, 'subgingival_landscape.pdf'));
end

function plot_tree_edges(principal_tree, edges)
[row_index, column_index] = find(triu(edges) > 0);
for i = 1:numel(row_index)
    node_pair = [row_index(i), column_index(i)];
    plot3(principal_tree(1,node_pair), principal_tree(2,node_pair), ...
        principal_tree(3,node_pair), '-k', 'LineWidth', 1.2);
end
end
