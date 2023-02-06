function [MS, CPL, GE, MCIC, MCOC, MBC, M, A] = getNetworkDirectedMeasures( ...
    connectivity_matrix, eegChannelsName, th)


%% measures obtained from the graphs ***************************************
% step 1: define a threshold to discard weak and non-significant connections
% Rubinov2010: "Threshold values are often arbitrarily determined, and
% networks should ideally be characterized across a broad range of
% thresholds."


n_chans = numel(eegChannelsName);


%% directed measures obtained with Matlab functions:

connectivity_matrix_th = connectivity_matrix;
connectivity_matrix_th(abs(connectivity_matrix)>=th);
G_th = digraph(connectivity_matrix, eegChannelsName, 'omitselfloops');


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

% mean incloseness centrality
MCIC = mean(centrality(G_th,'incloseness')); % from matlab

% mean outcloseness centrality
MCOC = mean(centrality(G_th,'outcloseness')); % from matlab

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

% MWCC = clustering_coef_wd(connectivity_matrix_th);

% prone to influences resulting from over-sampling and common sources


% transitivity: rate at which nodes with a common neighbor are connected

% T = transitivity_wu(connectivity_matrix_th);


% TO COMPUTE MWCC AND T IT IS REQUIRED THAT THE connectivity_matrix_th 
% CONTAINS all weights must be between 0 and 1

% classic modularity
gamma = 1;
if any(any(connectivity_matrix_th))
    [~, M] = modularity_dir(connectivity_matrix_th, gamma);
else
    M = 0;
end


%% Measures of resilience
% from Brain Connectivity Toolbox
% assortativity: tendency of nodes to connect to other nodes with similar
% properties

A = assortativity_wei(connectivity_matrix_th, 0);


end