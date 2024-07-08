# Financial Data Science project: Data Extraction for Results Presentation
# Computer Engineering: Intelligent Control Systems
# UniPV - 2023/24 
# student: Federico Pietro Bernardini
# ID: 514959

# This script is used to extract the authors of the papers for each year from the dataset.

import pandas as pd
import seaborn as sns

# Data reading
data = pd.read_csv('../dblp-v10.csv')

# Remove rows with NaN values
data = data.dropna()

# Just use the first 25000 rows for this example
data = data.head(25000)

# Filter data from 1980 to 2010
data = data[(data['year'] >= 1980) & (data['year'] <= 2009)]

# Just keep usefull columns
data = data[['authors', 'year']]

# Dividing data by year and then exctraction of the authors associated with each year
data_grouped = data.groupby('year')
dataframes = [(year, group) for year, group in data_grouped]

# For each dataframe group, create a list of authors for each year
authors_per_year = []
for year, df in dataframes:
    authors = [author for authors in df['authors'] for author in authors[1:-1].replace("'", "").split(", ")]
    authors_per_year.append((year, authors))

# save the results
authors_per_year = pd.DataFrame(authors_per_year, columns=['year', 'authors'])
authors_per_year.to_csv('../Datasets/authors.csv', index=False)