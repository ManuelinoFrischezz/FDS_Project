# Financial Data Science project: Data Cleaning & Preprocessing for Time Series Analysis
# Computer Engineering: Intelligent Control Systems
# UniPV - 2023/24 
# student: Federico Pietro Bernardini
# ID: 514959

# EXTRACTION OF THE ADJANCECY MATRICES FOR REFERNCES-PAPERS RELATIONSHIP

import pandas as pd

# Data reading
data = pd.read_csv('../dblp-v10.csv')

# Remove rows with NaN values
data = data.dropna()

# Just use the first 25000 rows for this example
data = data.head(25000)

# Filter data from 1980 to 2010
data = data[(data['year'] >= 1980) & (data['year'] <= 2009)]

# Keep only the necessary columns
data = data[['id', 'references', 'year']]

# Group data by the year of publication creating a dataframe for each year
data_grouped = data.groupby('year')

# Create a list of dataframes, one for each year
dataframes = [(year, group) for year, group in data_grouped]

# For each dataframe group, create an adjacency matrix for the references and ids
for year, df in dataframes:
    # Clean and split the 'references' column
    df['references'] = df['references'].apply(lambda x: x.strip('[]').replace("'", "").split(", ") if pd.notna(x) else [])

    # Create a sorted list of all unique ids and references
    id_list = df['id'].tolist()
    reference_list = [ref for refs in df['references'] for ref in refs]
    
    # Combine both lists to ensure all nodes are represented in the matrix
    all_papers = sorted(set(id_list + reference_list))

    # Creating the directed adjacency matrix with all_papers as both rows and columns
    adj_matrix = pd.DataFrame(0, index=all_papers, columns=all_papers)

    # Fill the adjacency matrix
    for idx, row in df.iterrows():
        paper_id = row['id']
        references = row['references']
        for ref in references:
            if ref in adj_matrix.columns:  # Ensure the reference is in the column list
                adj_matrix.loc[paper_id, ref] += 1

    # Save the adjacency matrix as a csv file
    adj_matrix.to_csv(f'../Datasets/TimeSeriesReferencePaper/adjacency_matrix_{year}.csv')
