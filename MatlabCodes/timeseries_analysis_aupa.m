%% Financial Data Science Project: Data Analysis of Bipartite Binary Networks
% Computer Engineering: Intelligent Control Systems
% UniPV - 2023/24 
% student: Federico Pietro Bernardini
% ID: 514959

% In this file the analysis over the Authors-Papers networks have been 
% carried out.

%% Cleaning the environment

clear
close all
clc

%% Import the Brain Connectivity

addpath('/Users/federicobernardini/Desktop/Project-FDS/2019_03_03_BCT')

%% Importing Networks Data (just bipartite adj matrix)

% use just one dataset at a time
folder = '../Datasets/TimeSeriesAuthorPaper/'; % folder with datasets
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

%% Author-Papers metrics evolution in year

degrees_aupa = cell(1,ndatasets);
r_aupa = zeros(1,ndatasets); % assortativity
dimensions_aupa = zeros(1,ndatasets);
avg_degrees_aupa = zeros(1,ndatasets);
betweenness_cen_aupa = cell(1,ndatasets); % betweenness centrality
mean_distance_aupa = zeros(1, ndatasets);

for i=1:ndatasets
    degrees_aupa{1,i} = degrees_und(A{i});
    avg_degrees_aupa(1,i) = mean(degrees_aupa{1,i}); 
    dimensions_aupa(1,i) = size(A{i},1)+size(A{i},2);
    r_aupa(1,i) = assortativity_bin(A_bipartite{i},0);
    betweenness_cen_aupa{1,i} = betweenness_bin(A_bipartite{i});

    distance = distance_bin(A_bipartite{i});
    distance(distance==inf) = nan; % substitute inf -> nan for, computing the mean nan are not considered
    no_nan_distance = distance(distance > 0 & ~isnan(distance)); % removing the NaN and diagonal elements (distance from themself) from the matrix
    mean_distance_aupa(1,i) = mean(no_nan_distance);
end
% r_aupa(isnan(r_aupa)) = 0; % nan values -> 0

save AuthorPaper/data

%% Single bipartite graph visualization

index = 10; % offset year of the graph we want to plot (1980 + offset year = actual year)

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

text = sprintf('Authors-Papers Graph (%i)', 1980+index-1);
title(text, 'Interpreter', 'latex', 'FontSize', 30)

hold on
s1 = scatter(nan, nan, 'r', 'filled', 'DisplayName', 'Authors'); % Red dot for first group
s2 = scatter(nan, nan, 'b', 'filled', 'DisplayName', 'Papers'); % Blue dot for second group
hold off
legend([s1, s2], 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);

%% Degree Distribution in 1980, 1995 & 2009 (AuthorsPapers)

load AuthorPaper/data.mat

studied_year = [1, 16, 30];
plot_elements = 1; % number of elements for the subplot

% Dimensions of the opened figure window
width = 800;  % width in pixels
height = 600; % height in pixels
left = 400;   % horizontal position (distance from the left corner of the screen)
bottom = 200; % vertical position (distance from the bottom corner of the screen)

figure('Position', [left, bottom, width, height]);
for i = 1:length(studied_year)
    
    % Power-law fitting for set U (Authors)
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

    % Creating the log-log plot for set U
    subplot(3, 2, plot_elements)
    scatter(bin_centers, normalized_degree, 'blue', 'filled', 'o', 'SizeData', 50); 
    grid on; hold on;
    plot(x, k(1) * x.^-k(2), "LineWidth", 2);

    ylabel('$P(k)$', 'Interpreter', 'latex', 'FontSize', 25)
    xlabel('$k$', 'Interpreter', 'latex', 'Fontsize', 25)
    text = sprintf('$\\mathcal{U}$ log-log Degree Distribution %d', years(studied_year(i)));
    title(text, 'Interpreter', 'latex', 'FontSize', 30)
    str = sprintf('$P(k) = ck^{-\\alpha}$, $c=%.2f, \\alpha=%.2f$', k(1), k(2));
    legend('Raw Data', str, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'southwest'); 
    xticks(sort(0:1:max(degreesU{1, studied_year(i)})))
    set(gca, 'Xscale', 'log', 'YScale', 'log');
    ylim([0 1])
    xtickangle(45)
    
    %%%%%%%%%%%%%%%%%
    
    % Poissonian fitting for set V (papers)
    x = -1:max(degreesV{1,studied_year(i)}); % x-axis values for the Poisson distribution
    poisson_pdf = poisspdf(x, avg_degrees_aupa(1,studied_year(i)));
    
    % Smoothing the poisson distribution
    x_interp = -1:0.1:max(degrees_aupa{1,studied_year(i)});
    poisson_pdf_interp = interp1(x, poisson_pdf, x_interp, 'spline');

    % Taking data from the histogram for the fitting
    h = histogram(degreesV{1, studied_year(i)}, 'Normalization', 'pdf', 'Visible', 'off', 'HandleVisibility', 'off');
    normalized_degree = h.Values;
    bin_edges = h.BinEdges;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

    % Estimation of the best parameters for set V
    initial_guess = [0 0];
    k = lsqcurvefit(model, initial_guess, bin_centers, normalized_degree);
    x = linspace(0.5, max(bin_centers), 100);

    % Creating the log-log plot
    subplot(3,2,plot_elements+1)
    scatter(bin_centers, normalized_degree, 'blue', 'filled', 'o', 'SizeData',50); 
    grid on; hold on;
    plot(x_interp, poisson_pdf_interp, 'red', 'LineWidth',2)
    
    ylabel('$P(k)$', 'Interpreter', 'latex', 'FontSize', 25)
    xlabel('$k$','Interpreter', 'latex', 'Fontsize', 25)
    text = sprintf('$\\mathcal{V}$ log-log Degree Distribution %d', years(studied_year(i)));
    title(text, 'Interpreter', 'latex', 'FontSize', 30)
    text = sprintf('$P(k) = \\frac{\\mu^k}{k!} e^{-\\mu}$, $\\mu = %.2f$', avg_degrees_aupa(1,studied_year(i)));
    legend('Raw Data',text, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'southwest')
    xticks(sort([1 0:5:max(degrees_aupa{1,studied_year(i)})]))
    set(gca, 'XScale', 'log', 'YScale', 'log');


    plot_elements = plot_elements + 2;

end

%% Top 5 authors in terms of papers published

load General/names_per_year.mat
load AuthorPaper/data.mat

raw_authors_per_years = cell(1,ndatasets);
authors_per_years = cell(1,ndatasets);

for i=1:ndatasets

    raw_authors_per_years{i} = names_per_year{1,i}(:,'authors');
    
    % Extraction the single names and removing useless symbols
    single_names = [];
    for j=1:size(raw_authors_per_years{i},1)
        text = table2cell(raw_authors_per_years{i}(j,1)); 
        splitted_names = split(text{1}(3:end-2), ',');
        single_names = [single_names splitted_names'];
    end
    
    % Final names per year structure
    authors_per_years{i} = single_names;
    
    % Needs to clean this variable every year considered
    clear single_names
end

% Counting the most influential authors per year
authors_count_per_year = cell(1, ndatasets);

for i = 1:ndatasets
    % Extract unique authors and count occurrences
    unique_authors = unique(authors_per_years{i});
    num_occurrences = zeros(1, length(unique_authors));
    
    for j = 1:length(unique_authors)
        num_occurrences(j) = sum(strcmp(authors_per_years{i}, unique_authors{j}));
    end
    
    % Create the new structure
    authors_count_per_year{i} = [unique_authors; num2cell(num_occurrences)];
    
    % Sort columns based on the number of occurrences (second row)
    [~, sort_idx] = sort(cell2mat(authors_count_per_year{i}(2,:)), 'descend');
    authors_count_per_year{i} = authors_count_per_year{i}(:, sort_idx);
    
    if size(authors_count_per_year{i}, 2) > 5
        authors_count_per_year{i} = authors_count_per_year{i}(:, 1:5);
    end
end

% Bar plot

% Create a map that links author to a color:
% Collect all author names through years to create a unique list
all_authors = cell(1,ndatasets*5);
aux_authors = []; % just an auxiliary variables
offset = 0;
for i = 1:ndatasets
    aux_authors = [aux_authors authors_count_per_year{i}(1,1:5)];
    for j=1:5
        all_authors{j+offset} = authors_count_per_year{i}(1,j);
    end
    offset = offset + 5;
end
all_unique_authors = unique(aux_authors, 'stable');

% Generate a colormap with unique colors
num_authors = length(all_unique_authors);
color_map = jet(num_authors);

% Create a map from author names to colors
author_color_map = dictionary(all_unique_authors, num2cell(color_map, 2)');

% Plotting the result
offset = 0;
space_between_groups = 1; % space between the groups of bars
group_width = 5; % number of bars per group

figure()
grid on
hold on
for i = 1:ndatasets
    x_positions = (1:group_width) + offset;

    b = bar(x_positions, [authors_count_per_year{1,i}{2,:}], 'FaceColor', 'flat');

    colors = zeros(5, 3);
    % Adding the color through a map
    for j = 1:5        
        author_name = all_authors{(i-1)*5 + j};
        colors(j, :) = cell2mat(author_color_map(author_name));
    end 
    b.CData = colors;

    offset = offset + group_width + space_between_groups;
end

xticks(3:5+space_between_groups:200)
xticklabels(years)
xlabel('year', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$k$', 'Interpreter', 'latex', 'FontSize', 25)
title('Top 5 Authors By Annual Publications', 'Interpreter', 'latex', 'FontSize', 30)

% Legend creation
legend_entries = gobjects(1, length(all_unique_authors));
hold on
for k = 1:length(all_unique_authors)
    legend_entries(k) = patch(NaN, NaN, cell2mat(author_color_map(all_unique_authors(k))));
end

legend(legend_entries, all_unique_authors, 'Location', 'southoutside', 'FontSize', 10, 'NumColumns', 6);