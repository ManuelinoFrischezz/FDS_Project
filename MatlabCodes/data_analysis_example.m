%% Financial Data Science Project: BC Toolbox understanding
% Computer Engineering: Intelligent Control Systems
% UniPV - 2023/24 
% student: Federico Pietro Bernardini
% ID: 514959

% ATTENTION: this code works only for binary undirected symmetric adj matrix.
% Its only porpouse is to undetrstand how the BCT works, no results from
% here has been used in the final report.

%% Cleaning the environment

clear
close all
clc

%% Import the Brain Connectivity

addpath('/Users/federicobernardini/Desktop/Project-FDS/2019_03_03_BCT')
addpath('/Users/federicobernardini/Desktop/Project-FDS/2019_03_03_BCT/data_and_demos')
addpath('/Users/federicobernardini/Desktop/Project-FDS/CountingNodesStructures')

%% Importing Network Data

filename = '../Datasets/example_adjacency_matrix.csv';
A = readmatrix(filename);
A = A(:,2:end); % adjacency matrix

%% Plotting the Graph

% Dimensions of the opened figure window
width = 800;  % width in pixels
height = 600; % hight in pixels
left = 400;   % horizontal position (distance from the left corner of the screen)
bottom = 200; % vertical position (distance from the bottom corner of the screen)

figure('Position', [left, bottom, width, height]);

G = graph(A);
degrees = degrees_und(A); % degrees of each node

Gplot = plot(G);
layout(Gplot, 'force')
Gplot.NodeCData = degrees;
nsizes = 2*sqrt(degrees-min(degrees)+0.2);
nsizes(nsizes < 3) = 3;
Gplot.MarkerSize = nsizes;
Gplot.LineWidth = 2;
% Gplot.MarkerSize = 2; % use a big value just if you want to see a sub-network
% highlight(Gplot,2419,'MarkerSize',20) % to highlight a specific node
c = colorbar;
colormap('jet')
ylabel(c, '$k$', 'FontSize', 20, 'Rotation', 0, 'Interpreter','latex');

title('Authors-Authors Network (1990)', 'Interpreter', 'latex', 'FontSize', 30)

%% Degree Distribution

mean_degree = mean(degrees); % lambda for the poissonian distribution

% Poisson pdf given my data
x = -1:max(degrees); % x-axis values for the Poisson distribution
poisson_pdf = poisspdf(x, mean_degree);

% Smoothing the poisson distribution
x_interp = -1:0.1:max(degrees);
poisson_pdf_interp = interp1(x, poisson_pdf, x_interp, 'spline');

% Real distribution vs estimated ones
figure('Position', [left, bottom, width, height]);
subplot(1,2,1)
histogram(degrees,'Normalization','pdf'); 
grid on; hold on
plot(x_interp, poisson_pdf_interp, 'r', 'LineWidth', 2); 
hold on
xline(mean_degree, '--', 'LineWidth', 2, 'Color', 'm')

ylabel('$P(k)$', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('$k$','Interpreter', 'latex', 'Fontsize', 25)
title('Degree Distribution', 'Interpreter', 'latex', 'FontSize', 30)
legend('Raw Data','Poissonian Distribution','Mean Value', 'Interpreter', 'latex', 'FontSize', 20)
xticks(0:2:max(degrees))
ylim([0 0.3])
% xtickangle(30)

% Creating the log-log plot
h = histogram(degrees, 'Normalization', 'pdf', 'Visible', 'off', 'HandleVisibility', 'off');
normalized_degree = h.Values;
bin_edges = h.BinEdges;
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;

subplot(1,2,2)
scatter(bin_centers, normalized_degree, 'blue', 'filled', 'o', 'SizeData',50); 
grid on; hold on;
plot(x_interp,poisson_pdf_interp, 'red', 'LineWidth',2)

ylabel('$P(k)$', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('$k$','Interpreter', 'latex', 'Fontsize', 25)
title('log-log Degree Distribution', 'Interpreter', 'latex', 'FontSize', 30)
legend('Raw Data','Poissonian Distribution', 'Interpreter', 'latex', 'FontSize', 20)
xticks(sort([1 0:5:max(degrees)]))
xlim([1 max(degrees)])
ylim([10e-5 1])
set(gca, 'XScale', 'log', 'YScale', 'log');

%% Paths-length

distance = distance_bin(A); % it takes some time, de-comment just if you need it

distance(distance==inf) = nan; % substitute inf -> nan for, computing the mean nan are not considered
no_nan_distance = distance(distance > 0 & ~isnan(distance)); % removing the NaN and diagonal elements (distance from themself) from the matrix
mean_distance = mean(no_nan_distance);

% Power-law fitting
model = @(b, x) x.^(-b(1)); % model = x^(-α) (Power Law)

% Taking data from the istogram for the fitting
figure()
h = histogram(no_nan_distance,'Normalization','pdf');
counts = h.Values;
bin_edges = h.BinEdges;
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;

% Estimation of the best parameters
initial_guess = 1;
k = lsqcurvefit(model, initial_guess, bin_centers, counts);
x = linspace(0.5, max(bin_centers), 100);

% Data distribution vs power law distribution
grid on; hold on
plot(x,x.^-k,"LineWidth", 2);
hold on
xline(mean_distance, '--', 'LineWidth', 2, 'Color', 'm')

ylabel('$P(d)$', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('$d$','Interpreter', 'latex', 'Fontsize', 25)
title('Distance Distribution', 'Interpreter', 'latex', 'FontSize', 30)
str = sprintf('$P(d) = d^{-\\alpha}$, $\\alpha=%.4f$', k);
legend('Raw Data', str, 'Mean Value', 'Interpreter', 'latex', 'FontSize', 20);
xticks(-1:max(max(distance)))
xlim([0.5 max(max(distance))])
ylim([0 1])
axis square

%% Cluestering Coefficient (NOT WORKING)

C = clustering_coef_bu(A); % clustering coefficient for each node

% Plot of the clustering coefficient distribution
figure()
h = histogram(C,'Normalization','probability','BinEdges',0:0.05:2.1); grid on;
bin_edges = h.BinEdges;
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;

ylabel('$P(C)$', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('$C$','Interpreter', 'latex', 'Fontsize', 25)
title('Clustering Coefficient Distribution', 'Interpreter', 'latex', 'FontSize', 30)
xticks(bin_centers)

% Plot of the C(k)-k loglog plot
% figure()
% scatter(degrees,C,'blue', 'filled', 'o', 'SizeData',50); grid on


% xticks(1:2:max(degrees))
% set(gca, 'XScale', 'log', 'YScale', 'log');
% axis square

%% Centralities

% Degree 
Dc = degrees/(size(A,1)-1);

% Betweenness
Bc = betweenness_bin(A);
Bc = Bc/((size(A,1)-1)*(size(A,1)-2)); % normalization

% Eigenvector
Ec = eigenvector_centrality_und(A);

figure()
Gplot = plot(G);
layout(Gplot, 'force')
Gplot.NodeCData = Dc;
Gplot.LineWidth = 2;
nsizes = 10*sqrt(Dc-min(Dc)+0.2);
Gplot.MarkerSize = nsizes;
% Gplot.MarkerSize = 10; % use a big value just if you want to see a sub-network
% highlight(Gplot,2419,'MarkerSize',20) % to highlight a specific node
title('Degree Centrality','Interpreter', 'latex', 'FontSize', 30)
colorbar

figure()
Gplot = plot(G);
layout(Gplot, 'force')
Gplot.NodeCData = Bc;
Gplot.LineWidth = 2;
nsizes = 10*sqrt(Bc-min(Bc)+0.2);
Gplot.MarkerSize = nsizes;
% Gplot.MarkerSize = 10; % use a big value just if you want to see a sub-network
% highlight(Gplot,2419,'MarkerSize',20) % to highlight a specific node
title('Betweenness Centrality','Interpreter', 'latex', 'FontSize', 30)
colorbar

figure()
Gplot = plot(G);
layout(Gplot, 'force')
Gplot.NodeCData = Ec;
Gplot.LineWidth = 2;
nsizes = 10*sqrt(Ec-min(Ec)+0.2);
Gplot.MarkerSize = nsizes;
% Gplot.MarkerSize = 10; % use a big value just if you want to see a sub-network
% highlight(Gplot,2419,'MarkerSize',20) % to highlight a specific node
title('Eigenvector Centrality','Interpreter', 'latex', 'FontSize', 30)
colorbar

%% Assortativity

r = assortativity_bin(A,0);

% e_jk = p_j * p_k (if no degree correlation)
E = repmat(normalized_degree,size(normalized_degree,2),1);
E = E.*E';

% Assortativity matrix 
figure()
imagesc(E);

colormap("pink");
xlabel('$k$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$j$', 'FontSize', 20, 'Interpreter', 'latex');
c = colorbar;
ylabel(c, '$e_{jk}$', 'FontSize', 20, 'Rotation', 0, 'Interpreter','latex');
title('$e_{jk} = p_j * p_k$','Interpreter', 'latex', 'FontSize', 30)
axis square

%% Average Next Neighbor Degree

annd = zeros(1,size(degrees,2));
for i=1:size(A,1)
    annd(1,i) = (A(i,:)*degrees')/degrees(i);
end
annd(isnan(annd)) = 0; % nan -> 0, beacause are the ones with k=0 (no neighbors)

% Average annd for each degree k
unique_degrees = unique(degrees); % different values of degrees 0,1,...,19
mean_annd = zeros(size(unique_degrees));

for i = 2:length(unique_degrees)
    % indeces for unique degree i
    indices = degrees == unique_degrees(i);
    % mean of the elements with index in indeces vector
    mean_annd(i) = mean(annd(indices));
end

% Model estimation
model = @(b, x) b(1)*x.^b(2); % model = ax^(γ)
initial_guess = [0 0];
params = lsqcurvefit(model, initial_guess, degrees, annd, counts);
x = linspace(0,max(degrees)+1,1000);

% Plot the results
figure()
scatter(degrees, annd, "blue", "filled", "o", "LineWidth", 2); grid on; hold on
plot(x, params(1)*x.^params(2), "m--", "LineWidth", 2); hold on
scatter(unique_degrees, mean_annd, "red", "filled", "o", "LineWidth", 3)

xlabel('$k$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$K_{annd}(k)$', 'FontSize', 20, 'Interpreter', 'latex');
title('Average Next Neighbor Degree','Interpreter', 'latex', 'FontSize', 30)
% str = sprintf('$K_{annd}(k)$ = %.4f$k^{%.4f}$', params(1), params(2));
str = sprintf('$K_{annd}(k)$ = $ak^{\\gamma}$ ($\\gamma=%.2f$)', params(2));
legend('Raw Data',str, '$<K_{annd}(k)>$', 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northwest')
xlim([1 20])
set(gca, 'XScale', 'log', 'YScale', 'log');
axis square

%% Communities

gamma = 1;
% [Ci,Q] = modularity_und(A, gamma); % less precise in communities identification
[Ci,Q] = community_louvain(A,gamma);

num_community = max(Ci);

% Assigning different colors to different communities
colors = lines(num_community);
node_colors = zeros(size(Ci,1), 3);

for i = 1:size(Ci,1)
    community_idx = Ci(i);
    node_colors(i, :) = colors(community_idx, :);
end

% Plot the results
figure()
h = plot(G, 'NodeColor', node_colors, 'MarkerSize', 7, 'NodeLabel', {});

title('Graph Communities','Interpreter', 'latex', 'FontSize', 30)

% Dendogram is useless, all the clusters are disconnected

%% Motifs: Degree-Preserving Random Graph

ITER = 2; % rewiring parameter, each edge is rewired approximately ITER times
% [A_rand,~] = randmio_und_connected(A,ITER); % decomment only if you really need it, really long computational time

% Plot the results
figure()
% subplot(1,2,1)
% Gplot = plot(G);
% Gplot.NodeCData = degrees;
% Gplot.LineWidth = 2;
% nsizes = 2*sqrt(degrees-min(degrees)+0.2);
% Gplot.MarkerSize = nsizes;
% title('Original Graph','Interpreter', 'latex', 'FontSize', 30)
% c = colorbar;
% ylabel(c, '$k$', 'FontSize', 20, 'Rotation', 0, 'Interpreter','latex');

% subplot(1,2,2)
degrees_rand = degrees_und(A_rand);
Gplot = plot(graph(A_rand));
Gplot.NodeCData = degrees_rand;
Gplot.LineWidth = 2;
nsizes = 2*sqrt(degrees_rand-min(degrees_rand)+0.2);
Gplot.MarkerSize = nsizes;
title('Degree-Preserving Randomisation Network','Interpreter', 'latex', 'FontSize', 30)
c = colorbar;
ylabel(c, '$k$', 'FontSize', 20, 'Rotation', 0, 'Interpreter','latex');

%% Motifs: 3 Nodes Structures

[~,F3] = motif3struct_bin(A);
[~, F3_rand] = motif3struct_bin(A_rand);

% Used motifs (only motifs 9 and 13 are used)
mstruct1 = find_motif34(9,3);
mstruct2 = find_motif34(13,3);

pdf3 = sum(F3,2)/sum(F3,"all");

% Motifs visualization
figure()
subplot(1,2,1)
m = graph(mstruct1(:,:,1));
mplot = plot(m);
layout(mplot, 'force')
mplot.MarkerSize = 10;
mplot.LineWidth = 2;
mplot.NodeColor = 'red';
mplot.NodeFontSize = 15;
title('Motif 9','Interpreter', 'latex', 'FontSize', 25)
str = sprintf('$P(m=9) =$ %.2f', pdf3(9));
legend(str, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'southeast')

subplot(1,2,2)
m = graph(mstruct2);
mplot = plot(m);
layout(mplot, 'force')
mplot.MarkerSize = 10;
mplot.LineWidth = 2;
mplot.NodeColor = 'red';
mplot.NodeFontSize = 15;
title('Motif 13','Interpreter', 'latex', 'FontSize', 25)
str = sprintf('$P(m=13) =$ %.2f', pdf3(13));
legend(str, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'southeast')

% Motifs distribution per node
x3 = linspace(1,size(A,1),size(A,1));
y3 = sum(F3,1);

figure()
bar(x3, y3, 10, 'grouped', 'blue'); hold on; grid on
ylabel('$\#$ of motifs', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('$n$','Interpreter', 'latex', 'Fontsize', 25)
title('3 Node-Motifs Distribution', 'Interpreter', 'latex', 'FontSize', 30)

% Z-score computation
triangles_per_node_rand = sum(F3_rand,1); % each element is the number of complete triangles that a node is part of
triangles_per_node = sum(F3,1); % tot number of triangles in the orginal graph

zscore3 = zeros(1,size(triangles_per_node,2));
for i=1:size(triangles_per_node,2)
    zscore3(i) = (triangles_per_node(i) - mean(triangles_per_node_rand))/std(triangles_per_node_rand);
end

% Z-score results
figure()
histogram(zscore3, 'Normalization','probability'); grid on; hold on;
xline(mean(zscore3), 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--')

ylabel('$P(Z)$', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('$Z$','Interpreter', 'latex', 'Fontsize', 25)
title('Z-score PDF (3 nodes motifs)', 'Interpreter', 'latex', 'FontSize', 30)
legend('', 'Mean Value', 'Interpreter', 'latex', 'FontSize', 20)
axis square

%% Motifs: 4 Nodes Structures

[~,F4] = motif4struct_bin(A);
[~, F4_rand] = motif4struct_bin(A_rand);

motifs = find(sum(F4,2));

% Used motifs (only motifs 9 and 13 are used)
mstruct1 = find_motif34(113,4);
mstruct2 = find_motif34(126,4);
mstruct3 = find_motif34(128,4);
mstruct4 = find_motif34(199,4);

pdf4 = sum(F4,2)/sum(F4,"all");

% Motifs visualization
figure()
subplot(2,2,1)
m = graph(mstruct1(:,:,1));
mplot = plot(m);
layout(mplot, 'force')
mplot.MarkerSize = 10;
mplot.LineWidth = 2;
mplot.NodeColor = 'red';
mplot.NodeFontSize = 15;
title('Motif 113','Interpreter', 'latex', 'FontSize', 25)
str = sprintf('$P(m=113) =$ %.2f', pdf4(113));
legend(str, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'southeast')

subplot(2,2,2)
m = graph(mstruct2(:,:,1));
mplot = plot(m);
layout(mplot, 'force')
mplot.MarkerSize = 10;
mplot.LineWidth = 2;
mplot.NodeColor = 'red';
mplot.NodeFontSize = 15;
title('Motif 126','Interpreter', 'latex', 'FontSize', 25)
str = sprintf('$P(m=126) =$ %.2f', pdf4(126));
legend(str, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeast')

subplot(2,2,3)
m = graph(mstruct3(:,:,1));
mplot = plot(m);
layout(mplot, 'force')
mplot.MarkerSize = 10;
mplot.LineWidth = 2;
mplot.NodeColor = 'red';
mplot.NodeFontSize = 15;
title('Motif 128','Interpreter', 'latex', 'FontSize', 25)
str = sprintf('$P(m=128) =$ %.2f', pdf4(128));
legend(str, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northwest')

subplot(2,2,4)
m = graph(mstruct4);
mplot = plot(m);
layout(mplot, 'force')
mplot.MarkerSize = 10;
mplot.LineWidth = 2;
mplot.NodeColor = 'red';
mplot.NodeFontSize = 15;
title('Motif 199','Interpreter', 'latex', 'FontSize', 25)
str = sprintf('$P(m=199) =$ %.2f', pdf4(199));
legend(str, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'southeast')

% Motifs distribution per node
x4 = linspace(1,size(A,1),size(A,1));
y4 = sum(F4,1);

figure()
bar(x4, y4, 10, 'grouped', 'blue'); hold on; grid on
ylabel('$\#$ of motifs', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('$n$','Interpreter', 'latex', 'Fontsize', 25)
title('4 Node-Motifs Distribution', 'Interpreter', 'latex', 'FontSize', 30)

% Z-score computation
squares_per_node_rand = sum(F4_rand,1); % each element is the number of "complete square" that a node is part of
squares_per_node = sum(F4,1); % tot number of "complete square" in the orginal graph

zscore4 = zeros(1,size(squares_per_node,2));
for i=1:size(squares_per_node,2)
    zscore4(i) = (squares_per_node(i) - mean(squares_per_node_rand))/std(squares_per_node_rand);
end

% Z-score results
figure()
histogram(zscore4, 'Normalization','probability'); grid on;
xline(mean(zscore4), 'LineWidth', 2, +'Color', 'r', 'LineStyle', '--')

ylabel('$P(Z)$', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('$Z$','Interpreter', 'latex', 'Fontsize', 25)
title('Z-score PDF (4 nodes motifs)', 'Interpreter', 'latex', 'FontSize', 30)
legend('', 'Mean Value', 'Interpreter', 'latex', 'FontSize', 20)
axis square