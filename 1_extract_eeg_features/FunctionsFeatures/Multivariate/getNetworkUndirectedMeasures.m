function [MD, MS, CPL, GE, MCC, MBC, WGCC, T, M, A] = getNetworkUndirectedMeasures( ...
    connectivity_matrix, connectivity_vector, eegChannelsName, electrode_pairs, th)



% References
% (1) Lehnertz2017
% (2) Rubinov2010
% (3) Cohen2014


%% measures obtained from the graphs ***************************************
% step 1: define a threshold to discard weak and non-significant connections
% Rubinov2010: "Threshold values are often arbitrarily determined, and
% networks should ideally be characterized across a broad range of
% thresholds."


n_chans = numel(eegChannelsName);


%% undirected measures obtained with Matlab functions:

connections_no_node_name = [electrode_pairs, connectivity_vector];
connections_threshold = connections_no_node_name(connectivity_vector<th,:);
G_th = graph(connections_threshold(:,1), connections_threshold(:,2), ...
    connections_threshold(:,3), n_chans);


connectivity_matrix_th = connectivity_matrix;
connectivity_matrix_th(connectivity_matrix<th) = 0;

% set the diagonal to zero (Rubinov2010):
connectivity_matrix_th(find(eye(size(connectivity_matrix_th)))) = 0;
connect_mat_bin = double(connectivity_matrix_th~=0);


% average degree
% degree of an individual node is equal to the number of links connected to
% that node
% The mean network degree is most commonly used as a measure of density, or
% the total "wiring cost" of the network.
% the weighted variant of the degree, sometimes termed the strength, is
% defined as the sum of all neighboring link weights.
MD = mean(sum(connect_mat_bin)); % slightly faster method
% D = mean(degree(G_th)) % other method

% strength
MS = mean(sum(connectivity_matrix_th));


%% Measures of integration
% from Brain Connectivity Toolbox

% convert a weighted connection matrix to a weighted connection-length matrix
L = weight_conversion(connectivity_matrix, 'lengths');
% get the lengths of shortest paths between all pairs of nodes.
D = distance_wei(L);
% characteristic path length
[CPL, efficiency, ~, ~, ~] = charpath(D);

% global efficiency
GE = efficiency_wei(connectivity_matrix_th);

%% Measures of centrality

% mean closeness centrality
% closeness centrality: average of the distances (i.e., shortest paths) of
% a node to all other nodes
MCC = mean(centrality(G_th,'closeness')); % from matlab



% mean betweenness centrality
% betweenness centrality: ratio of all shortest paths (between all nodes)
% passing through a given node.

% BC = betweenness_wei(L);
MBC = mean(centrality(G_th, 'betweenness')); % from matlab


%% Measures of segregation
% from Brain Connectivity Toolbox

% mean weighted clustering coefficient
% clustering coefficient a node is the ratio of that node's neighbors that
% are connected to each other (known as the local clustering coefficient)
% get the weighted global clustering coefficient
WGCC = mean(clustering_coef_wu(connectivity_matrix_th));
% prone to influences resulting from over-sampling and common sources


% transitivity: rate at which nodes with a common neighbor are connected
T = transitivity_wu(connectivity_matrix_th);


% classic modularity
gamma = 1;
if any(any(connectivity_matrix_th))
    [~, M] = modularity_und(connectivity_matrix_th, gamma);
else
    M = 0;
end

% clustering coefficient entropy (NO LITERATURE FOR THIS MEASURE besides Abbas2021)
% C = clustering_coef_bu(connect_mat_bin)
%
% [c] = unique(C);
%
% prob = zeros(numel(c),1);
% for ii = 1:length(c)
%     prob(ii) = sum(C==c(ii));
% end
% nonz_prob = prob(find(prob));
% CCE = -sum(nonz_prob.*log2(nonz_prob));


%% Measures of resilience
% from Brain Connectivity Toolbox
% assortativity: tendency of nodes to connect to other nodes with similar
% properties

A = assortativity_wei(connectivity_matrix_th, 0);


%% small-worldness
% This measure may be useful for snapshot characterization of an ensemble
% of networks, but it may also falsely report a small- world topology in
% highly segregated, but poorly integrated networks. Consequently, this
% measure should not in general be regarded as a substitute for individual
% assessments of integration and segregation. Rubinov2010


end