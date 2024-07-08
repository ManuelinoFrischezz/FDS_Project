%% Financial Data Science Project: Data Analysis of Bipartite Binary Networks
% Computer Engineering: Intelligent Control Systems
% UniPV - 2023/24 
% student: Federico Pietro Bernardini
% ID: 514959

% In this file the analysis over the Journals-Papers networks have been 
% carried out.

%% Cleaning the environment

clear
close all
clc

%% Import the Brain Connectivity

addpath('/Users/federicobernardini/Desktop/Project-FDS/2019_03_03_BCT')

%% Importing Networks Data (just bipartite adj matrix)

% use just one dataset at a time
folder = '../Datasets/TimeSeriesPaperJournal/'; % folder with datasets
filelist = dir(fullfile(folder, '*.csv')); % file in the folder

% Years extraction
years = zeros(1,length(filelist));
for i=1:length(filelist)
    years(i) = str2double(regexp([filelist(i).name], '\d+', 'match'));
end

% Ordering the files by the year
[sortedYears, sortIdx] = sort(years);
sortedFilelist = filelist(sortIdx);

ndatasets = numel(filelist);% n of used datasets
A = cell(1,ndatasets); % Structure with all the datasets
for i = 1:ndatasets
    filename = fullfile(folder, filelist(i).name); % file complete name
    A{i} = readmatrix(filename);
    A{i} = A{i}(:,2:end); % adjacency matrix
end

ndatasets = length(years);

%% Network Structure Creation

G_bipartite = cell(1,ndatasets);
A_bipartite = cell(1,ndatasets);
degreesU = cell(1,ndatasets);
degreesV = cell(1,ndatasets);
labels = cell(1,ndatasets);

for i=1:ndatasets

    % Remove nodes with no connections
    degreesU{1,i} = sum(A{i}, 2); % Degree of nodes in U
    degreesV{1,i} = sum(A{i}, 1); % Degree of nodes in V

    % Find indices of nodes with connections
    nodesU = find(degreesU{1,i} > 0);
    nodesV = find(degreesV{1,i} > 0);

    % Create new adjacency matrix without isolated nodes
    A{i} = A{i}(nodesU, nodesV);

    % N° of nodes of the two sets
    nU = size(A{i},1); % first set size of the bipartite network
    nV = size(A{i},2); % second set size of the bipartite network
    
    % Nodes labels
    % labelsU = arrayfun(@(x) ['U', num2str(x)], 1:nU, 'UniformOutput', false);
    % labelsV = arrayfun(@(x) ['V', num2str(x)], 1:nV, 'UniformOutput', false);
    % labels{i} = [labelsU, labelsV];
    
    % Adding the nodes of the second set as additional nodes in the adj matrices
    A_bipartite{i} = [zeros(nU, nU), A{i}; 
                   A{i}', zeros(nV, nV)];
    
    % Creating the bipartite graph
    G_bipartite{i} = graph(A_bipartite{i});

end

%% Paper-Journal metrics evolution in year

degrees_pajo = cell(1,ndatasets);
r_pajo = zeros(1,ndatasets); % assortativity
dimensions_pajo = zeros(1,ndatasets);
avg_degrees_pajo = zeros(1,ndatasets);
betweenness_cen_pajo = cell(1,ndatasets); % betweenness centrality
mean_distance_pajo = zeros(1,ndatasets);

for i=1:ndatasets
    avg_degrees_pajo(1,i) = mean(degreesU{i}); 
    dimensions_pajo(1,i) = size(A{i},1)+size(A{i},2);
    r_pajo(1,i) = assortativity_bin(A_bipartite{i},0);
    betweenness_cen_pajo{1,i} = betweenness_bin(A_bipartite{i});

    distance = distance_bin(A_bipartite{i});
    distance(distance==inf) = nan; % substitute inf -> nan for, computing the mean nan are not considered
    no_nan_distance = distance(distance > 0 & ~isnan(distance)); % removing the NaN and diagonal elements (distance from themself) from the matrix
    mean_distance_pajo(1,i) = mean(no_nan_distance);
end

save PaperJournal/data 

%% Load data

load PaperJournal/data 

%% Single bipartite graph visualization

index = 26; % offset year of the graph we want to plot (1980 + offset year = actual year)

% N° of nodes of the two sets
nU = size(A{index},1); % first set size of the bipartite network
nV = size(A{index},2); % second set size of the bipartite network

% Dimensions of the opened figure window
width = 800;  % width in pixels
height = 600; % hight in pixels
left = 400;   % horizontal position (distance from the left corner of the screen)
bottom = 200; % vertical position (distance from the bottom corner of the screen)

figure('Position', [left, bottom, width, height]);

Gplot = plot(G_bipartite{index});
% NOTE: reduce the second MarkerSize to better display Papers-Journal Graph
% for big graphs
highlight(Gplot, 1:nU, 'NodeColor', 'r','MarkerSize', 5);
highlight(Gplot, nU+1:nU+nV, 'NodeColor', 'b', 'MarkerSize', 5); 
% labelnode(Gplot, 1:nU+nV, labels{index});

text = sprintf('Papers-Journals Graph (%i)', 1980+index-1);
title(text, 'Interpreter', 'latex', 'FontSize', 30)

% Adding a legend for Papers-Journals Graph
hold on
s1 = scatter(nan, nan, 'r', 'filled', 'DisplayName', 'Journals'); % Red dot for first group
s2 = scatter(nan, nan, 'b', 'filled', 'DisplayName', 'Papers'); % Blue dot for second group
hold off
legend([s1, s2], 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);

%% Degree Distribution in 1980, 1995 & 2009 (PaperJournals)

load PaperJournals/data.mat

% for set U (= Journals) in graph PaperJournals, because each paper is
% published just on one journal, so degree is always 1

studied_year = [1, 16, 30];
plot_elements = 1; % number of elements for the subplot

% Dimensions of the opened figure window
width = 800;  % width in pixels
height = 600; % height in pixels
left = 400;   % horizontal position (distance from the left corner of the screen)
bottom = 200; % vertical position (distance from the bottom corner of the screen)

figure('Position', [left, bottom, width, height]);
for i = 1:length(studied_year)
    
    % Power-law fitting
    model = @(b, x) b(1) * x.^(-b(2)); % model = cx^(-α) (Power Law)

    % Taking data from the histogram for the fitting
    h = histogram(degreesU{1, studied_year(i)}, 'Normalization', 'pdf', 'Visible', 'off', 'HandleVisibility', 'off');
    normalized_degree = h.Values;
    bin_edges = h.BinEdges;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

    % Estimation of the best parameters for set U
    initial_guess = [0 0];
    k = lsqcurvefit(model, initial_guess, bin_centers, normalized_degree);
    x = linspace(0.5, max(bin_centers), 100);
    
    % Creating the plot for set U
    subplot(3, 2, plot_elements)
    histogram(degreesU{1, studied_year(i)}, 'Normalization', 'pdf'); hold on;
    scatter(bin_centers, normalized_degree, 'blue', 'filled', 'o', 'SizeData', 50); 
    grid on; hold on;
    plot(x, k(1) * x.^-k(2), "LineWidth", 2, 'Color', 'red');

    ylabel('$P(k)$', 'Interpreter', 'latex', 'FontSize', 25)
    xlabel('$k$', 'Interpreter', 'latex', 'Fontsize', 25)
    text = sprintf('$\\mathcal{U}$ Degree Distribution %d', years(studied_year(i)));
    title(text, 'Interpreter', 'latex', 'FontSize', 30)
    str = sprintf('$P(k) = ck^{-\\alpha}$, $c=%.2f, \\alpha=%.2f$', k(1), k(2));
    legend('Raw Data', '', str, 'Interpreter', 'latex', 'FontSize', 20); 
    xticks(sort([1 0:10:max(degreesU{1, studied_year(i)})]))
    ylim([0 1])
    xtickangle(45)

    % Creating the log-log plot for set U
    subplot(3, 2, plot_elements+1)
    scatter(bin_centers, normalized_degree, 'blue', 'filled', 'o', 'SizeData', 50); 
    grid on; hold on;
    plot(x, k(1) * x.^-k(2), "LineWidth", 2);

    ylabel('$P(k)$', 'Interpreter', 'latex', 'FontSize', 25)
    xlabel('$k$', 'Interpreter', 'latex', 'Fontsize', 25)
    text = sprintf('$\\mathcal{U}$ log-log Degree Distribution %d', years(studied_year(i)));
    title(text, 'Interpreter', 'latex', 'FontSize', 30)
    str = sprintf('$P(k) = ck^{-\\alpha}$, $c=%.2f, \\alpha=%.2f$', k(1), k(2));
    legend('Raw Data', str, 'Interpreter', 'latex', 'FontSize', 20); 
    xticks(sort([1 0:10:max(degreesU{1, studied_year(i)})]))
    set(gca, 'Xscale', 'log', 'YScale', 'log');
    ylim([0 1])
    xtickangle(45)

    plot_elements = plot_elements + 2;

end

%% Top 5 journals in terms of papers published on them

load General/names_per_year.mat
load PaperJournal/data.mat

journals_per_year = cell(1,ndatasets);

for i=1:ndatasets
    journals_per_year{i} = names_per_year{1,i}{:,3};
end

% Counting the most influential authors per year
journals_count_per_year = cell(1, ndatasets);

for i = 1:ndatasets
    % Extract unique journals and count occurrences
    unique_journals = unique(journals_per_year{i}, 'stable');
    num_occurrences = zeros(1, length(unique_journals));
    
    for j = 1:length(unique_journals)
        num_occurrences(j) = sum(strcmp(journals_per_year{i}, unique_journals{j}));
    end
    
    % Create the new structure
    journals_count_per_year{i} = [unique_journals'; num2cell(num_occurrences)];
    
    % Sort columns based on the number of occurrences (second row)
    [~, sort_idx] = sort(cell2mat(journals_count_per_year{i}(2,:)), 'descend');
    journals_count_per_year{i} = journals_count_per_year{i}(:, sort_idx);
    
    if size(journals_count_per_year{i}, 2) > 5
        journals_count_per_year{i} = journals_count_per_year{i}(:, 1:5);
    end
end

% Bar plot

% Create a map that links author to a color:
% Collect all author names through years to create a unique list
all_journals = cell(1,ndatasets*5);
aux_journals = []; % just an auxiliary variables
offset = 0;
for i = 1:ndatasets
    aux_journals = [aux_journals journals_count_per_year{i}(1,1:5)];
    for j=1:5
        all_journals{j+offset} = journals_count_per_year{i}(1,j);
    end
    offset = offset + 5;
end
all_unique_journals = unique(aux_journals, 'stable');

% Generate a colormap with unique colors
num_journals = length(all_unique_journals);
color_map = jet(num_journals);

% Create a map from author names to colors
journals_color_map = dictionary(all_unique_journals, num2cell(color_map, 2)');

% Plotting the result
offset = 0;
space_between_groups = 1; % space between the groups of bars
group_width = 5; % number of bars per group

figure()
grid on
hold on
for i = 1:ndatasets
    x_positions = (1:group_width) + offset;

    b = bar(x_positions, [journals_count_per_year{1,i}{2,:}], 'FaceColor', 'flat');

    colors = zeros(5, 3);
    % Adding the color through a map
    for j = 1:5        
        author_name = all_journals{(i-1)*5 + j};
        colors(j, :) = cell2mat(journals_color_map(author_name));
    end 
    b.CData = colors;

    offset = offset + group_width + space_between_groups;
end

xticks(3:5+space_between_groups:200)
xticklabels(years)
xlabel('year', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$k$', 'Interpreter', 'latex', 'FontSize', 25)
title('Top 5 Journal By Annual Publications', 'Interpreter', 'latex', 'FontSize', 30)

% Legend creation
legend_entries = gobjects(1, length(all_unique_journals));
hold on
for k = 1:length(all_unique_journals)
    legend_entries(k) = patch(NaN, NaN, cell2mat(journals_color_map(all_unique_journals(k))));
end

legend(legend_entries, all_unique_journals, 'Location', 'southoutside', 'FontSize', 13, 'NumColumns', 2);