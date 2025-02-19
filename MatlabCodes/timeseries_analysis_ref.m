%% Financial Data Science Project: Data Analysis of Directed Binary Networks
% Computer Engineering: Intelligent Control Systems
% UniPV - 2023/24 
% student: Federico Pietro Bernardini
% ID: 514959

% In this file the analysis over the References-Papers network have been 
% carried out.

%% Cleaning the environment

clear
close all
clc

%% Import the Brain Connectivity

addpath('/Users/federicobernardini/Desktop/Project-FDS/2019_03_03_BCT')

%% Importing Networks Data (just directed symmetric adj matrix)

% ATTENTION: JUST LOAD THE DATA (see later)

folder = '../Datasets/TimeSeriesReferencePaper/'; % folder with datasets
filelist = dir(fullfile(folder, '*.csv')); % file in the folder

% Years extraction
years = zeros(1, length(filelist));
for i = 1:length(filelist)
    years(i) = str2double(regexp([filelist(i).name], '\d+', 'match'));
end

% Ordering the files by the year
[sortedYears, sortIdx] = sort(years);
sortedFilelist = filelist(sortIdx);

ndatasets = numel(filelist); % n of used datasets
A = cell(1, ndatasets); % Structure with all the datasets
for i = 1:ndatasets
    filename = fullfile(folder, filelist(i).name); % file complete name
    A{i} = readmatrix(filename);
    A{i} = A{i}(2:end, 2:end); % adjacency matrix (use it for Reference-Paper Graph)

    % Remove isolated nodes
    degrees = sum(A{i}, 1) + sum(A{i}, 2)'; % sum of rows and columns
    nodes = find(degrees > 0); % nodes with at least one connection
    A{i} = A{i}(nodes, nodes); % Keep only nodes with connections
end

ndatasets = length(years);

% save 'ReferencePaperData/ReferencePaperData' -v7.3 % too long computational time, better save the data and load them

%% Loading Data

load ReferencePaperData/ReferencePaperData.mat

%% Network Structure Creation & Analysis

G = cell(1,ndatasets);
degrees = cell(1,ndatasets);
avg_in_degrees = zeros(1,ndatasets); % avg in degree
avg_out_degrees = zeros(1,ndatasets); % avg out degree
avg_degrees_ref = zeros(1,ndatasets); % avg degree
dimensions_ref = zeros(1,ndatasets); % number of nodes
r_ref = zeros(4,ndatasets); % assortativity
betweenness_cen = cell(1,ndatasets); % betweenness centrality
mean_distance_ref = zeros(1, ndatasets);

for i=1:ndatasets
    G{i} = digraph(A{i});
    [id,od,deg] = degrees_dir(A{i}); % degrees of each node (id = inner degree, od = outer degree, deg = id+od)
    degrees{1,i} = id; % id = inner degree
    degrees{2,i} = od; % od = outer degree
    degrees{3,i} = deg; % deg = id+od
    avg_in_degrees(1,i) = mean(degrees{1,i});
    avg_out_degrees(1,i) = mean(degrees{2,i});
    avg_degrees_ref(1,i) = mean(degrees{3,i});
    dimensions_ref(1,i) = size(A{i},1);
    r_ref(1,i) = assortativity_bin(A{i},1); % out-degree/in-degree correlation
    r_ref(2,i) = assortativity_bin(A{i},2); % in-degree/out-degree correlation
    r_ref(3,i) = assortativity_bin(A{i},3); % out-degree/out-degree correlation
    r_ref(4,i) = assortativity_bin(A{i},4); % in-degree/in-degree correlation
    betweenness_cen{1,i} = betweenness_bin(A{i});
    betweenness_cen{1,i} = betweenness_cen{1,i}/((size(A,1)-1)*(size(A,1)-2)); % to normalize the values

    distance = distance_bin(A{i});
    distance(distance==inf) = nan; % substitute inf -> nan for, computing the mean nan are not considered
    no_nan_distance = distance(distance > 0 & ~isnan(distance)); % removing the NaN and diagonal elements (distance from themself) from the matrix
    mean_distance_ref(1,i) = mean(no_nan_distance);

end
r_ref(isnan(r_ref)) = 0; % nan values -> 0

% save ReferencePaperData/data

%% Load data

load ReferencePaperData/data.mat

%% Single directed graph visualization

index = 11; % offset year of the graph we want to plot (1980 + offset year - 1 = actual year)

% Dimensions of the opened figure window
width = 800;  % width in pixels
height = 600; % hight in pixels
left = 400;   % horizontal position (distance from the left corner of the screen)
bottom = 200; % vertical position (distance from the bottom corner of the screen)

figure('Position', [left, bottom, width, height]);

Gplot = plot(G{index});
Gplot.NodeCData = degrees{3,index};
nsizes = 3*sqrt(degrees{3,index}-min(degrees{3,index})+0.2);
% nsizes(nsizes < 3) = 3;
Gplot.MarkerSize = nsizes;
Gplot.LineWidth = 2;
Gplot.ArrowSize = 10;
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
height = 600; % height in pixels
left = 400;   % horizontal position (distance from the left corner of the screen)
bottom = 200; % vertical position (distance from the bottom corner of the screen)

figure('Position', [left, bottom, width, height]);
for i = 1:length(studied_year)
    
    % Power-law fitting
    model = @(b, x)b(1) * x.^(-b(2)); % model = cx^(-α) (Power Law)
    
    % Taking data from the histogram for the fitting
    h = histogram(degrees{3, studied_year(i)}, 'Normalization', 'pdf', 'Visible', 'off', 'HandleVisibility', 'off');
    counts = h.Values;
    bin_edges = h.BinEdges;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    
    % Estimation of the best parameters
    initial_guess = [0 0];
    k = lsqcurvefit(model, initial_guess, bin_centers, counts);
    x = linspace(0.5, max(bin_centers), 100);
    
    % Creating the normal plot
    h = histogram(degrees{3, studied_year(i)}, 'Normalization', 'pdf', 'Visible', 'off', 'HandleVisibility', 'off');
    normalized_degree = h.Values;
    bin_edges = h.BinEdges;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    
    subplot(3, 2, plot_elements)
    histogram(degrees{3, studied_year(i)}, 'Normalization', 'pdf');grid on
    hold on
    scatter(bin_centers, normalized_degree, 'blue', 'filled', 'o', 'SizeData', 50)
    hold on
    plot(x, k(1) * x.^-k(2), "LineWidth", 2, "Color", 'red');
    
    ylabel('$P(k)$', 'Interpreter', 'latex', 'FontSize', 25)
    xlabel('$k$', 'Interpreter', 'latex', 'Fontsize', 25)
    text = sprintf('Degree Distribution %d', years(studied_year(i)));
    title(text, 'Interpreter', 'latex', 'FontSize', 30)
    str = sprintf('$P(k) = ck^{-\\alpha}$, $c=%.2f, \\alpha=%.2f$', k(1), k(2));
    legend('','Raw Data', str, 'Interpreter', 'latex', 'FontSize', 20); 
    xticks(sort([1 0:10:max(degrees{3, studied_year(i)})]))
    ylim([0 1])
    xtickangle(45)
   

    % Creating the log-log plot
    h = histogram(degrees{3, studied_year(i)}, 'Normalization', 'pdf', 'Visible', 'off', 'HandleVisibility', 'off');
    normalized_degree = h.Values;
    bin_edges = h.BinEdges;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    
    subplot(3, 2, plot_elements+1)
    scatter(bin_centers, normalized_degree, 'blue', 'filled', 'o', 'SizeData', 50); 
    grid on; hold on;
    plot(x, k(1) * x.^-k(2), "LineWidth", 2);
    
    ylabel('$P(k)$', 'Interpreter', 'latex', 'FontSize', 25)
    xlabel('$k$', 'Interpreter', 'latex', 'Fontsize', 25)
    text = sprintf('log-log Degree Distribution %d', years(studied_year(i)));
    title(text, 'Interpreter', 'latex', 'FontSize', 30)
    str = sprintf('$P(k) = ck^{-\\alpha}$, $c=%.2f, \\alpha=%.2f$', k(1), k(2));
    legend('Raw Data', str, 'Interpreter', 'latex', 'FontSize', 20); 
    xticks(sort([1 0:10:max(degrees{3, studied_year(i)})]))
    set(gca, 'Xscale', 'log', 'YScale', 'log');
    ylim([0 1])
    xtickangle(45)

    plot_elements = plot_elements + 2;

end

%% Assortativity

% Dimensions of the opened figure window
width = 800;  % width in pixels
height = 600; % hight in pixels
left = 400;   % horizontal position (distance from the left corner of the screen)
bottom = 200; % vertical position (distance from the bottom corner of the screen)

f = figure('Position', [left, bottom, width, height]);
plot(years, r_ref(1,:), "LineWidth", 2); grid on
hold on
plot(years, r_ref(2,:), "LineWidth", 2)
hold on
plot(years, r_ref(3,:), "LineWidth", 2)
hold on
plot(years, r_ref(4,:), "LineWidth", 2)

ylabel('$r^\rightarrow$', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('$year$', 'Interpreter', 'latex', 'Fontsize', 25)
title('References-Papers Assortativity', 'Interpreter', 'latex', 'FontSize', 30)
legend('out-in degree', 'in-out degree', 'out-out degree', 'in-in degree', 'Interpreter', 'latex', 'FontSize', 20); 
xticks(sort(0:2:max(years)))
xlim([min(years) max(years)])
ylim([-1.1 1.1])
xtickangle(45)

%% Betweenness Centrality

n = 30; % analyzed year

% Find top 5 referenced papers
positions = find(betweenness_cen_{n});
positions = sort(positions, "descend");

if length(positions) > 5
    positions = positions(1:5);
end

%
figure()
Gplot = plot(G{n});
layout(Gplot, 'force')
Gplot.NodeCData = betweenness_cen{n};
Gplot.LineWidth = 2;
nsizes = sqrt(betweenness_cen{n} - min(betweenness_cen{n})+0.2);
Gplot.MarkerSize = nsizes;
% Gplot.MarkerSize = 10; % use a big value just if you want to see a sub-network
% highlight(Gplot, positions, 'MarkerSize', 20) % to highlight a specific node
title('Betweenness Centrality','Interpreter', 'latex', 'FontSize', 30)
colorbar

%% Top 5 in & out degrees publications
% (Just in degrees results are displayed)

load ReferencePaperData/data.mat
load General/names_per_year

% Finding the names of the papers divided by year

% inner degree = more cited papers
% outer degree = more citing papers

names = cell(4,ndatasets); 
% 1st row = in-degree value
% 2nd row = in-degree id
% 3rd row = out-degree value
% 4th row = out-degree id

for i=1:ndatasets

    % in-degree top 5 selection
    id = [degrees{1,i}(degrees{1,i} > 0); find(degrees{1,i} > 0)];

    % Sorting the elements
    [~, sort_idx] = sort(id(1,:), 'descend');
    sorted_id = id(:, sort_idx);

    names{1,i} = sorted_id(1,:); % values
    names{2,i} = sorted_id(2,:); % id
    if length(names{1,i}) > 3
        names{1,i} = names{1,i}(1:5);
        names{2,i} = names{2,i}(1:5);
    end

    % out-degree top 5 selection
    id = [degrees{2,i}(degrees{2,i} > 0); find(degrees{2,i} > 0)];
   
    % Sorting the elements
    [~, sort_idx] = sort(id(1,:), 'descend');
    sorted_id = id(:, sort_idx);

    names{3,i} = sorted_id(1,:); % values
    names{4,i} = sorted_id(2,:); % id
    if length(names{3,i}) > 5 
        names{3,i} = names{3,i}(1:5);
        names{4,i} = names{4,i}(1:5);
    end

end

titles = cell(2,ndatasets); % row 1 = indegree; row 2 = outdegree;
prova = names;

for i=1:ndatasets
    
    prova{2,i}(prova{2,i} > size(names_per_year{1,i},1)) = nan;
    prova{2,i}(isnan(prova{2,i})) = 0;
    prova{4,i}(prova{4,i} > size(names_per_year{1,i},1)) = nan;
    prova{4,i}(isnan(prova{2,i})) = 0;
    
    for j=1:length(prova{2,i})
        if ~isnan(prova{2,i}(j)) & prova{2,i}(j) > 0
            titles{1,i}{j} = names_per_year{1,i}(prova{2,i}(j), "title"); % in degree papers name
        else
            titles{1,i}{j} = cell2table({'unknown'});
        end

        if ~isnan(prova{4,i}(j)) & prova{4,i}(j) > 0
            titles{2,i}{j} = names_per_year{1,i}(prova{4,i}(j), "title"); % out degree papers name
        else
            titles{2,i}{j} = cell2table({'unknown'});
        end

    end
end

% Creating the table with all the titles
papers_names = [];
for i=1:size(titles,2)
    for j=1:5
        if istable(titles{1,i}{j})
            papers_names = [papers_names table2cell(titles{1,i}{j})];
        else
            papers_names = [papers_names 'unknown'];
        end
    end
end

% Creation of the list for papers title with the right order of occurence
ordered_paper_titles = cell(1,length(unique(papers_names)));

% Initializing the cell randomly
for i=1:length(ordered_paper_titles)
    ordered_paper_titles{i} = 'mario'; % random name just to be sure is not the name of a paper
end

j = 1;
for i = 1:length(papers_names)
    if ~ismember(papers_names{i}, ordered_paper_titles)
        ordered_paper_titles{j} = papers_names{i};
        j = j+1;
    end
end

% Generating a random colours scale for bars
num_colors = length(ordered_paper_titles);
colors = jet(num_colors);

% Creating a map to associate titles to colors
paper_colors = containers.Map();
for i = 1:length(ordered_paper_titles)
    if strcmp(ordered_paper_titles{i}, 'unknown')
        paper_colors(ordered_paper_titles{i}) = [0, 0, 0]; % black color -> 'unknown'
    else
        paper_colors(ordered_paper_titles{i}) = colors(i, :);
        % if two papers have the same name then have the same color 
    end
end

% Ectracting the unique paper names and colors from the map
unique_paper_titles = keys(paper_colors);
unique_paper_colors = values(paper_colors);

% Plotting the result
offset = 0;
space_between_groups = 1; % space between the groups of bars
group_width = 5; % number of bars per group

figure()
grid on
hold on
for i = 1:ndatasets

    x_positions = (1:group_width) + offset;

    b = bar(x_positions, names{1,i}, 'FaceColor', 'flat');
    
    colors = zeros(5, 3);
    % Assigning the colors through the map values
    for j = 1:5
        paper_title = papers_names{(i-1)*5 + j};
        colors(j, :) = paper_colors(paper_title);
    end
    b.CData = colors;

    offset = offset + group_width + space_between_groups;
end

xticks(3:5+space_between_groups:200)
xticklabels(years)
xlabel('year', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$k^{in}$', 'Interpreter', 'latex', 'FontSize', 25)
title('Top 5 Cited Papers By Year', 'Interpreter', 'latex', 'FontSize', 30)

% Legend creation
legend_entries = gobjects(1, length(ordered_paper_titles));
hold on
for k = 1:length(ordered_paper_titles)
    legend_entries(k) = patch(NaN, NaN, paper_colors(ordered_paper_titles{k}));
end

legend(legend_entries, ordered_paper_titles, 'Location', 'southoutside', 'FontSize', 10, 'NumColumns', 2);
