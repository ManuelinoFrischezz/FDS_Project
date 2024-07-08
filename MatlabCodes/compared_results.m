%% Financial Data Science Project: Timeseries Data Comparison 
% Computer Engineering: Intelligent Control Systems
% UniPV - 2023/24 
% student: Federico Pietro Bernardini
% ID: 514959

% In this file the plots to compare the graphs have been performed.
%% Cleaning the environment

clear
close all
clc

%% Degree evolution in years

load AuthorAuthorData/data.mat
load ReferencePaperData/data.mat
load PaperJournal/data.mat
load AuthorPaper/data.mat

% Dimensions of the opened figure window
width = 800;  % width in pixels
height = 600; % hight in pixels
left = 400;   % horizontal position (distance from the left corner of the screen)
bottom = 200; % vertical position (distance from the bottom corner of the screen)

%% Plot of the average degree time evolution

figure('Position', [left, bottom, width, height]);
plot(years, avg_degrees_aut, 'LineWidth', 2); grid on
hold on
plot(years, avg_degrees_ref, 'LineWidth', 2)
hold on
plot(years, avg_degrees_pajo, 'LineWidth', 2)
hold on
plot(years, avg_degrees_aupa, 'LineWidth', 2)

ylabel('$<k>$', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('year','Interpreter', 'latex', 'Fontsize', 25)
title('Average Degree Evolution', 'Interpreter', 'latex', 'FontSize', 30)
legend('Authors-Authors', 'References-Papers', 'Papers-Journals', 'Authors-Papers', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);
xlim([min(years) max(years)])
xticks(min(years):2:max(years))
xtickangle(30)

%% Plot of the nodes number time evolution

f = figure('Position', [left, bottom, width, height]);
plot(years, dimensions_aut, 'LineWidth', 2); grid on
hold on
plot(years, dimensions_ref, 'LineWidth', 2)
hold on
plot(years, dimensions_pajo, 'LineWidth', 2)
hold on
plot(years, dimensions_aupa, 'LineWidth', 2)

ylabel('$N$', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('year','Interpreter', 'latex', 'Fontsize', 25)
title('Nodes Evolution', 'Interpreter', 'latex', 'FontSize', 30)
legend('Authors-Authors', 'References-Papers', 'Papers-Journals', 'Authors-Papers', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);
xlim([min(years) max(years)])
xticks(min(years):max(years))
xtickangle(30)

%% Plot of the assortativity time evolution

figure('Position', [left, bottom, width, height]);
plot(years, r_aut, 'LineWidth', 2); grid on;
hold on
plot(years, r_pajo, 'LineWidth', 2)
hold on
plot(years, r_aupa, 'LineWidth', 2)
hold on
yline(0, 'LineWidth', 2, 'LineStyle','--', 'Color', 'red')

ylabel('$r$', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('year','Interpreter', 'latex', 'Fontsize', 25)
title('Assortativity Time Evolution', 'Interpreter', 'latex', 'FontSize', 30)
legend('Authors-Authors', 'Papers-Journals', 'Authors-Papers', '', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);
ylim([-1 1])
xlim([min(years) max(years)])
xticks(min(years):max(years))
xtickangle(30)

%% Plot the average path length time evolution 

% Dimensions of the opened figure window
width = 800;  % width in pixels
height = 600; % hight in pixels
left = 400;   % horizontal position (distance from the left corner of the screen)
bottom = 200; % vertical position (distance from the bottom corner of the screen)

figure('Position', [left, bottom, width, height]);
plot(years, mean_distance_aut, 'LineWidth', 2); grid on
hold on
plot(years, mean_distance_ref, 'LineWidth', 2)
hold on
plot(years, mean_distance_pajo, 'LineWidth', 2)
hold on
plot(years, mean_distance_aupa, 'LineWidth', 2)

ylabel('$<d>$', 'Interpreter', 'latex', 'FontSize', 25)
xlabel('year','Interpreter', 'latex', 'Fontsize', 25)
title('Average Path Length Evolution', 'Interpreter', 'latex', 'FontSize', 30)
legend('Authors-Authors', 'References-Papers', 'Papers-Journals', 'Authors-Papers', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);
xlim([min(years) max(years)])
xticks(min(years):2:max(years))
xtickangle(30)