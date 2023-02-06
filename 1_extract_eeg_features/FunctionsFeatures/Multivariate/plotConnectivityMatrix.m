function [] = plotConnectivityMatrix(connectivity_matrix, ...
    connectivity_vector, eegChannelsName, electrode_pairs, th, ...
    name_connectivity)


n_chans = numel(eegChannelsName);

figure()
subplot(121)
% connectivity_matrix(find(eye(size(connectivity_matrix)))) = 0;
imagesc(connectivity_matrix)
xticks(1:n_chans)
yticks(1:n_chans)
xticklabels(eegChannelsName)
yticklabels(eegChannelsName)
colorbar
xtickangle(90)
set(gca,'xaxisLocation','top')
title(regexprep(name_connectivity,'_',' '))

subplot(122)
if isempty(connectivity_vector)
    % directed connectivity
    if any(any(connectivity_matrix<0))
        connectivity_matrix(abs(connectivity_matrix)<th) = eps;
    else
        connectivity_matrix(connectivity_matrix<th) = eps;
    end
    G = digraph(connectivity_matrix, eegChannelsName, 'omitselfloops');
else
    % undirected connectivity
    connections2plot = [eegChannelsName(electrode_pairs(:,1))', ...
        eegChannelsName(electrode_pairs(:,2))' num2cell(connectivity_vector)];
    connections2plot(connectivity_vector<th,3) = num2cell(eps);
    
    G = graph(connections2plot(:,1), connections2plot(:,2), ...
        cell2mat(connections2plot(:,3)));
end

h = plot(G);
g_nodes = table2cell(G.Nodes);

positionsx = [4 6 3.8 6.2 3.5 6.5 3.8 6.2 4 6 3 7 2.5 7.5 3 7 5 5 5];

positionsy = [9.5 9.5 7.3 7.3 5 5 2.7 2.7 0.5 0.5 7.7 7.7 5 5 2.3 2.3 7.2 5 2.8];

x_changed = zeros(n_chans,1);
y_changed = x_changed;
for nn = 1:numel(g_nodes)
    ind_variable = strcmp(g_nodes(nn), eegChannelsName);
    x_changed(nn) = positionsx(ind_variable);
    y_changed(nn) = positionsy(ind_variable);
end

h = plot(G,'XData', x_changed, 'YData', y_changed);

weight_edges = G.Edges.Weight;
end_nodes = G.Edges.EndNodes;
weight_edges_unique = unique(weight_edges);

n_weights = numel(weight_edges_unique);
map = flipud(brewermap(n_weights,'Spectral'));
if any(weight_edges_unique==eps)
    map(weight_edges_unique==eps,:) = [0.75 0.75 0.75];
end

for ii = 1:n_weights
    indexes = weight_edges==weight_edges_unique(ii);
    nodes = end_nodes(indexes,:);
    for jj = 1:size(nodes,1)
        highlight(h,nodes(jj,:),'EdgeColor',map(ii,:))
    end
end

colormap(gca, map);

set(gca,'ytick',[])
set(gca,'xtick',[])

axis tight
% title(name_feature)

cb = colorbar;
ylabel(cb, regexprep(name_connectivity,'_',' ')) % , 'Interpreter', 'Latex', 'Fontsize', 12

if any(weight_edges_unique<0) % if there are negative values
    am = get(gca,'Alphamap');
    n_ticks = 20;
    cb.Ticks = linspace(min(am), max(am), n_ticks);
    downsampled = interp1(1:length(weight_edges_unique), weight_edges_unique, ...
        linspace(1,length(weight_edges_unique),n_ticks));
    cb.TickLabels = num2cell(round(downsampled*10)/10);
else
    if any(weight_edges_unique==eps)
        if ~all(weight_edges_unique==eps)
            % if there are non significant connections
            am = get(gca,'Alphamap');
            n_ticks = min([20 numel(weight_edges_unique)]);
            cb.Ticks = linspace(min(am), max(am), n_ticks);
            downsampled = interp1(1:length(weight_edges_unique), weight_edges_unique, ...
                linspace(1,length(weight_edges_unique),n_ticks));
            cb.TickLabels = num2cell(round(downsampled*10)/10);
        end
    else
        set(gca, 'clim', [weight_edges_unique(2) max(weight_edges_unique)]);
    end
end



if ~all(weight_edges_unique==eps)
    G.Edges.LWidths = 7*weight_edges/max(abs(weight_edges));
    h.LineWidth = abs(G.Edges.LWidths);
end
end