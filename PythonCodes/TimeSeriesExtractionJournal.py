# Financial Data Science project: Data Cleaning & Preprocessing for Time Series Analysis
# Computer Engineering: Intelligent Control Systems
# UniPV - 2023/24 
# student: Federico Pietro Bernardini
# ID: 514959

# EXTRACTION OF THE ADJANCECY MATRICES FOR PAPERS-JOURNALS RELATIONSHIP

import pandas as pd

# Data reading
data = pd.read_csv('../dblp-v10.csv')

# Remove rows with NaN values
data = data.dropna()

data = data.head(25000)

# Just data from 1980 to 2010 are used
data = data[(data['year'] >= 1980) & (data['year'] <= 2009)]

# Only keep the names of the authors, the title of the paper, the journal, and the year of publication.
data = data[['venue', 'id', 'year']]

# Group data by the year of publication creating a dataframe for each year
data_grouped = data.groupby('year')

# Create a list of dataframes, one for each year
dataframes = [(year, group) for year, group in data_grouped]

# For each dataframe group create an adjacency matrix for the venues and paper ids
for year, df in dataframes:
    # Create a sorted list of the venues and paper ids
    venues_list = df['venue'].dropna().unique().tolist()  # Remove NaN values from venues list
    ids_list = df['id'].unique().tolist()

    # Sort the lists
    venues_list = sorted(venues_list)
    ids_list = sorted(ids_list)

    # Create the bipartite adjacency matrix with venues as rows and paper ids as columns
    adj_matrix = pd.DataFrame(0, index=venues_list, columns=ids_list)

    # Fill the adjacency matrix
    for idx, row in df.iterrows():
        venue = row['venue']
        paper_id = row['id']
        if pd.isna(venue):
            continue  # Skip if venue is NaN
        adj_matrix.loc[venue, paper_id] = 1

    # Save the adjacency matrix as a csv file
    adj_matrix.to_csv(f'../Datasets/TimeSeriesPaperJournal/adjacency_matrix_{year}.csv')
