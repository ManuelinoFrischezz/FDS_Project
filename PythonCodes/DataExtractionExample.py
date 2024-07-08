# Financial Data Science project: Data Cleaning & Preprocessing
# Computer Engineering: Intelligent Control Systems
# UniPV - 2023/24 
# student: Federico Pietro Bernardini
# ID: 514959

# This script is used to preprocess the dataset and create the adjacency matrix for the authors of the papers.
# The dataset created is the one used just for understanding how to use the BCT, it has not been used for the analysis.

import pandas as pd

# data reading
data = pd.read_csv('../dblp-v10.csv')

n = 2000 # number of rows to use for the analysis

# I just use the first 5000 rows of the dataset
data = data.head(n)

# I just use data between 2000 and 2010
data = data[data['year'] >= 1980]
data = data[data['year'] <= 2009]

# print(data_n.head())

# just the names of the authors, the title of the paper, the number of citazions and the pubblication year are useful for the analysis
# removing column 1,4,5 and 8.
data = data.drop(data.columns[[0,3,4,7]], axis=1)

# create a sorted list of the authors
authors_list = []
for authors in data['authors']:
    authors = authors[1:-1].replace("'","").split(", ")
    for author in authors:
        authors_list.append(author)

# remove duplicates
authors_list = sorted(list(set(authors_list)))

# creating the adjacency matrix for the authors
adj_matrix = pd.DataFrame(0, index=authors_list, columns=authors_list)

# fill the adjacency matrix
for authors in data['authors']:
    authors = authors[1:-1].replace("'","").split(", ")
    for i in range(len(authors)):
        for j in range(i+1, len(authors)):
            adj_matrix[authors[i]][authors[j]] += 1
            adj_matrix[authors[j]][authors[i]] += 1

# save the adjacency matrix as a csv file
adj_matrix.to_csv('../example_adjacency_matrix.csv')

# print(adj_matrix.head())