# Financial Data Science project: Data Cleaning & Preprocessing for Time Series Analysis
# Computer Engineering: Intelligent Control Systems
# UniPV - 2023/24 
# student: Federico Pietro Bernardini
# ID: 514959

# EXTRACTION OF THE ADJANCECY MATRICES FOR AUTHORS-AUTHORS RELATIONSHIP

import pandas as pd

# Data reading
data = pd.read_csv('../dblp-v10.csv')

# Remove rows with NaN values
data = data.dropna()

data = data.head(25000)

# Just data from 1980 to 2009 are used
data = data[data['year'] >= 1980]
data = data[data['year'] <= 2009]

# Only keep the names of the authors and the title of the paper, and the year of publication.
# Assuming the columns are in the order: 'id', 'authors', 'title', 'journal', 'year', etc.
data = data[['authors', 'title', 'year']]

# Group data by the year of publication creating a dataframe for each year
data_grouped = data.groupby('year')

# Create a list of dataframes, one for each year
dataframes = []
for year, group in data_grouped:
    dataframes.append((year, group))

# For each dataframe group create an adjacency matrix for the authors co-authorship
for year, df in dataframes:
    # Create a sorted list of the authors
    authors_list = []
    for authors in df['authors']:
        # Assuming the authors are stored as a single string with authors separated by commas
        authors = authors[1:-1].replace("'", "").split(", ")
        for author in authors:
            authors_list.append(author)

    # Remove duplicates and sort the list
    authors_list = sorted(list(set(authors_list)))

    # Creating the adjacency matrix for the authors
    adj_matrix = pd.DataFrame(0, index=authors_list, columns=authors_list)

    # Fill the adjacency matrix
    for authors in df['authors']:
        authors = authors[1:-1].replace("'", "").split(", ")
        for i in range(len(authors)):
            for j in range(i + 1, len(authors)):
                adj_matrix.loc[authors[i], authors[j]] += 1
                adj_matrix.loc[authors[j], authors[i]] += 1

    # Save the adjacency matrix as a csv file
    adj_matrix.to_csv(f'../Datasets/TimeSeriesAuthorAuthor/adjacency_matrix_{year}.csv')