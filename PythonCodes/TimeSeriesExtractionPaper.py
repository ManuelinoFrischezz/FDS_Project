# Financial Data Science project: Data Cleaning & Preprocessing for Time Series Analysis
# Computer Engineering: Intelligent Control Systems
# UniPV - 2023/24 
# student: Federico Pietro Bernardini
# ID: 514959

# EXTRACTION OF THE ADJANCECY MATRICES FOR AUTHORS-PAPERS RELATIONSHIP

import pandas as pd

# Data reading
data = pd.read_csv('../dblp-v10.csv')

# Remove rows with NaN values
data = data.dropna()

data = data.head(25000)

# Just data from 1980 to 2009 are used
data = data[data['year'] >= 1980]
data = data[data['year'] <= 2009]

# Only keep the names of the authors, the title of the paper, the paper id, and the year of publication.
data = data[['authors', 'title', 'year', 'id']]

# Group data by the year of publication creating a dataframe for each year
data_grouped = data.groupby('year')

# Create a list of dataframes, one for each year
dataframes = []
for year, group in data_grouped:
    dataframes.append((year, group))

# For each dataframe group create an adjacency matrix for the authors and papers
for year, df in dataframes:
    # Create a sorted list of the authors and papers
    authors_list = []
    ids_list = df['id'].dropna().unique().tolist()  # Remove NaN values from journals list
    for authors in df['authors']:
        # Assuming the authors are stored as a single string with authors separated by commas
        authors = authors[1:-1].replace("'", "").split(", ")
        for author in authors:
            authors_list.append(author)

    # Remove duplicates and sort the list of authors
    authors_list = sorted(list(set(authors_list)))
    ids_list = sorted(ids_list)

    # Create the bipartite adjacency matrix with authors as rows and ids as columns
    adj_matrix = pd.DataFrame(0, index=authors_list, columns=ids_list)

    # Fill the adjacency matrix
    for idx, row in df.iterrows():
        authors = row['authors'][1:-1].replace("'", "").split(", ")
        id = row['id']
        if pd.isna(id):
            continue  # Skip if id is NaN
        for author in authors:
            adj_matrix.loc[author, id] += 1

    # Save the adjacency matrix as a csv file
    adj_matrix.to_csv(f'../Datasets/TimeSeriesAuthorPaper/adjacency_matrix_{year}.csv')