%% Financial Data Science Project: Convertion Table Building
% Computer Engineering: Intelligent Control Systems
% UniPV - 2023/24 
% student: Federico Pietro Bernardini
% ID: 514959

% In this file the creation of a table to convert the results from the
% graph to the real world meaning is done.

% Example of usage for Authors-Authors network: 
% Thanks to this table it is possible to understand which is the name of 
% the author that corresponds to node 28 in year 1998.

%% Cleaning the environment

clear
clc
close all

%% Creating the table of convertion to have the name of papers, journals and authors into the plots

ndatasets = 30;

convertion = readtable('../Datasets/convertion_table.csv');

% Estrai tutti gli anni unici nella tabella e ordinali
unique_years = unique(convertion.year);

% Inizializza una cell array per contenere i dati raggruppati per anno
names_per_year = cell(1, ndatasets);

% Riempi la cell array con i dati per ogni anno
for i = 1:ndatasets
    % Anno corrente
    current_year = unique_years(i);
    
    % Seleziona i dati per l'anno corrente
    data_for_year = convertion(convertion.year == current_year, :);
    
    % Memorizza i dati relativi all'anno corrente nella cell array
    names_per_year{i} = data_for_year;
end

save General/names_per_year names_per_year