%% Financial Data Science Project: Data Analysis of Undirected Binary Networks
% Computer Engineering: Intelligent Control Systems
% UniPV - 2023/24 
% student: Federico Pietro Bernardini
% ID: 514959

% In this file the analysis over the Authors-Authors network have been
% carried out.

%% Cleaning the environment

clear
close all
clc

%% Import the Brain Connectivity

addpath('/Users/federicobernardini/Desktop/Project-FDS/2019_03_03_BCT')

%% Importing Networks Data (just symmetric adj matrix)

folder = '../Datasets/TimeSeriesAuthorAuthor/'; % folder with datasets
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
    A{i} = A{i}(:,2:end); % adjacency matrix (use it for Author-Author Graph)
end

ndatasets = length(years);

%% Network Structure Creation

G = cell(1,ndatasets);
degrees_aut = cell(1,ndatasets);
avg_degrees_aut = zeros(1,ndatasets); % avg degree
dimensions_aut = zeros(1,ndatasets); % number of nodes
r_aut = zeros(1,ndatasets); % assortativity
betweenness_aut = cell(1,ndatasets);
mean_distance_aut = zeros(1,ndatasets);

for i=1:ndatasets
    G{i} = graph(A{i});
    degrees_aut{1,i} = degrees_und(A{i}); % degrees of each node
    avg_degrees_aut(1,i) = mean(degrees_aut{1,i}); 
    dimensions_aut(1,i) = size(A{i},1);
    r_aut(1,i) = assortativity_bin(A{i},0);
    betweenness_aut{1,i} = betweenness_bin(A{i});
    
    distance = distance_bin(A{i}); % it takes some time, de-comment just if you need it
    distance(distance==inf) = nan; % substitute inf -> nan for, computing the mean nan are not considered
    no_nan_distance = distance(distance > 0 & ~isnan(distance)); % removing the NaN and diagonal elements (distance from themself) from the matrix
    mean_distance_aut(1,i) = mean(no_nan_distance);
end
r_aut(isnan(r_aut)) = 0; % nan values -> 0

save AuthorAuthorData/data

%% Loading data

load AuthorAuthorData/data.mat

%% Single Graph Visualization

index = 11;

% Dimensions of the opened figure window
width = 800;  % width in pixels
height = 600; % hight in pixels
left = 400;   % horizontal position (distance from the left corner of the screen)
bottom = 200; % vertical position (distance from the bottom corner of the screen)

figure('Position', [left, bottom, width, height]);

Gplot = plot(G{index});
layout(Gplot, 'force')
Gplot.NodeCData = degrees_aut{index};
nsizes = 2*sqrt(degrees_aut{index}-min(degrees_aut{index})+0.2);
nsizes(nsizes < 3) = 3;
Gplot.MarkerSize = nsizes;
Gplot.LineWidth = 2;
% Gplot.MarkerSize = 2; % use a big value just if you want to see a sub-network
% highlight(Gplot,2419,'MarkerSize',20) % to highlight a specific node
c = colorbar;
colormap('jet')
ylabel(c, '$k$', 'FontSize', 20, 'Rotation', 0, 'Interpreter','latex');

text = sprintf('References-Papers Graph (%i)', 1980+index-1);
title(text, 'Interpreter', 'latex', 'FontSize', 30)

%% Degree Distribution in 1980, 1995 & 2009

studied_year = [1, 16, 30];
plot_elements = 1; % number of elements for the subplot

% Dimensions of the opened figure window
width = 800;  % width in pixels
height = 600; % hight in pixels
left = 400;   % horizontal position (distance from the left corner of the screen)
bottom = 200; % vertical position (distance from the bottom corner of the screen)

figure('Position', [left, bottom, width, height]);
for i=1:length(studied_year)
    
    % Poisson pdf given my data
    x = -1:max(degrees_aut{1,studied_year(i)})+1; % x-axis values for the Poisson distribution
    poisson_pdf = poisspdf(x, avg_degrees_aut(1,studied_year(i)));
    
    % Smoothing the poisson distribution
    x_interp = -1:0.1:max(degrees_aut{1,studied_year(i)})+1;
    poisson_pdf_interp = interp1(x, poisson_pdf, x_interp, 'spline');
    
    % Real distribution vs estimated ones
    subplot(3,2,plot_elements)
    h = histogram(degrees_aut{1,studied_year(i)}, 'Normalization', 'pdf');
    normalized_degree = h.Values;
    bin_edges = h.BinEdges;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;
    grid on; hold on
    plot(x_interp, poisson_pdf_interp, 'r', 'LineWidth', 2); 
    hold on
    xline(avg_degrees_aut(1,studied_year(i)), '--', 'LineWidth', 2, 'Color', 'm')
    
    ylabel('$P(k)$', 'Interpreter', 'latex', 'FontSize', 25)
    xlabel('$k$','Interpreter', 'latex', 'Fontsize', 25)
    text = sprintf('Degree Distribution %d', years(studied_year(i)));
    title(text, 'Interpreter', 'latex', 'FontSize', 30)
    text = sprintf('$\\langle k \\rangle = \\mu = %.2f$', avg_degrees_aut(1,studied_year(i)));
    legend('Raw Data', '$P(k) = \frac{\mu^k}{k!} e^{-\mu}$', text, 'Interpreter', 'latex', 'FontSize', 20)
    xticks(0:2:max(degrees_aut{1,studied_year(i)}))
    
    % Creating the log-log plot
    subplot(3,2,plot_elements+1)
    scatter(bin_centers, normalized_degree, 'blue', 'filled', 'o', 'SizeData',50); 
    grid on; hold on;
    plot(x_interp, poisson_pdf_interp, 'red', 'LineWidth',2)
    
    ylabel('$P(k)$', 'Interpreter', 'latex', 'FontSize', 25)
    xlabel('$k$','Interpreter', 'latex', 'Fontsize', 25)
    text = sprintf('log-log Degree Distribution %d', years(studied_year(i)));
    title(text, 'Interpreter', 'latex', 'FontSize', 30)
    text = sprintf('$P(k) = \\frac{\\mu^k}{k!} e^{-\\mu}$, $\\mu = %.2f$', avg_degrees_aut(1,studied_year(i)));
    legend('Raw Data',text, 'Interpreter', 'latex', 'FontSize', 20)
    xticks(sort([1 0:5:max(degrees_aut{1,studied_year(i)})]))
    set(gca, 'XScale', 'log', 'YScale', 'log');

    plot_elements = plot_elements + 2;

end

%% Most influential authors in terms of connections

load AuthorAuthorData/data.mat
load General/names_per_year.mat

% Find the name of the authors per year
authors = readtable('../Datasets/authors.csv');

authors_per_year = cell(1, size(authors, 1));
for i = 1:size(authors, 1)
    authors_string = authors{i, 2};
    authors_per_year{i} = split(authors_string{1}(2:end-2), ",")';
end

top5betweenness_per_year = cell(1, ndatasets);
indeces_per_year = cell(1, ndatasets); % number of nodes with betweenness > 0

% Sorting and selecting the top 5 values
for i=1:ndatasets
    indeces_per_year{i} = find(betweenness_aut{1,i});
    
    % aux structure: 1st row = values; 2nd row = nodes number; 3rd row = node-author name correspondance
    aux = [betweenness_aut{1,i}(1,indeces_per_year{i}); indeces_per_year{i}];
    [~, sortIdx] = sort(aux(2, :), 'descend'); % sorting aux in descend order
    aux = aux(:, sortIdx); % sorted aux
    
    % Selecting just the top 5 nodes
    if size(aux,2) > 5
        aux = aux(:,1:5);
    else
        % filling missing spaces with zeros
        for j=size(aux,2):4
            aux = [aux [0;0]];
        end
    end

    top5betweenness_per_year{i} = aux; 
end

% Selection of the 150 authors in the plot

for i=1:ndatasets
    for j=1:5
        if top5betweenness_per_year{1,i}(2,j) ~= 0
            all_authors{(i-1)*5+j} = authors_per_year{1,i}{1,top5betweenness_per_year{1,i}(2,j)};
        else
            all_authors{(i-1)*5+j} = 'unknown';
        end
    end
end
all_authors_unique = unique(all_authors, 'stable');

% Generate a colormap with unique colors
num_authors = length(all_authors_unique);
color_map = [[0 0 0]; jet(num_authors-1)]; % unknown element is black

% Create a map from author names to colors
author_color_map = dictionary(all_authors_unique, num2cell(color_map, 2)');

% Plot the result

offset = 0;
space_between_groups = 1; % space between the groups of bars
group_width = 5; % number of bars per group

figure()
grid on
hold on
for i = 1:ndatasets
    x_positions = (1:group_width) + offset;

    b = bar(x_positions, [top5betweenness_per_year{1,i}(1,:)], 'FaceColor', 'flat');

    colors = zeros(5, 3);
    % Adding the color through a map
    for j = 1:5        
        author_name = all_authors((i-1)*5 + j);
        colors(j, :) = cell2mat(author_color_map(author_name));
    end 
    b.CData = colors;


    offset = offset + group_width + space_between_groups;
end

xticks(3:5+space_between_groups:200)
xticklabels(years)
xlabel('year', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$B_C(\mathrm{author})$', 'Interpreter', 'latex', 'FontSize', 25)
title('Top 5 Annual Connecting Authors By Year', 'Interpreter', 'latex', 'FontSize', 30)

% Legend creation
legend_entries = gobjects(1, length(all_authors_unique));
hold on
for k = 1:length(all_authors_unique)
    legend_entries(k) = patch(NaN, NaN, cell2mat(author_color_map(all_authors_unique(k))));
end

legend(legend_entries, all_authors_unique, 'Location', 'southoutside', 'FontSize', 13, 'NumColumns', 4);