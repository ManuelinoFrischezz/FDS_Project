# Financial Data Science project: General Dataset Analysis
# Computer Engineering: Intelligent Control Systems
# UniPV - 2023/24 
# student: Federico Pietro Bernardini
# ID: 514959

# This script is used to perform a general analysis of the dataset and study it from an overall perspective.
# The analysis is divided into three parts:
# 1. Number of papers per year
# 2. Number of unique authors per year
# 3. Number of unique venues per year

import pandas as pd
import matplotlib.pyplot as plt

# Data reading
data = pd.read_csv('../dblp-v10.csv')

# Keeping only the necessary columns
data = data[['year', 'authors', 'venue']]

# Removing rows where 'authors' or 'venue' are NaN
data = data.dropna(subset=['authors', 'venue'])

# Counting the number of elements (papers) for each year
year_counts = data['year'].value_counts().sort_index()

# Function to clean and split the authors string
def clean_authors(authors):
    return authors[1:-1].replace("'", "").split(", ")

# Cleaning and splitting the 'authors' column
data['authors'] = data['authors'].apply(clean_authors)

# Exploding the 'authors' column to have one author per row
data_exploded = data.explode('authors')

# Counting the number of unique authors for each year
authors_per_year = data_exploded.groupby('year')['authors'].nunique()

# Counting the number of unique venues for each year
venues_per_year = data.groupby('year')['venue'].nunique()

# Plotting the bar chart for number of elements (papers) per year
plt.figure(figsize=(15, 7))
plt.bar(year_counts.index, year_counts.values)
plt.xticks(range(year_counts.index.min(), year_counts.index.max() + 1, 5))
plt.grid(True)
plt.axvline(1980, color="red", linewidth=2, linestyle="dashed")
plt.axvline(2009, color="red", linewidth=2, linestyle="dashed")
plt.xlabel('Year', fontsize=14)
plt.ylabel('Number of Papers', fontsize=14)
plt.title('Number of Papers per Year', fontsize=16)
plt.savefig('../Results/General/papers_per_year.png', dpi=300, bbox_inches='tight')
plt.show()

# Plotting the bar chart for number of unique authors per year
plt.figure(figsize=(15, 7))
plt.bar(authors_per_year.index, authors_per_year.values)
plt.xticks(range(authors_per_year.index.min(), authors_per_year.index.max() + 1, 5))
plt.grid(True)
plt.axvline(1980, color="red", linewidth=2, linestyle="dashed")
plt.axvline(2009, color="red", linewidth=2, linestyle="dashed")
plt.xlabel('Year', fontsize=14)
plt.ylabel('Number of Unique Authors', fontsize=14)
plt.title('Number of Unique Authors per Year', fontsize=16)
plt.savefig('../Results/General/authors_per_year.png', dpi=300, bbox_inches='tight')
plt.show()

# Plotting the bar chart for number of unique venues per year
plt.figure(figsize=(15, 7))
plt.bar(venues_per_year.index, venues_per_year.values)
plt.xticks(range(venues_per_year.index.min(), venues_per_year.index.max() + 1, 5))
plt.grid(True)
plt.axvline(1980, color="red", linewidth=2, linestyle="dashed")
plt.axvline(2009, color="red", linewidth=2, linestyle="dashed")
plt.xlabel('Year', fontsize=14)
plt.ylabel('Number of Unique Venues', fontsize=14)
plt.title('Number of Unique Venues per Year', fontsize=16)
plt.savefig('../Results/General/venues_per_year.png', dpi=300, bbox_inches='tight')
plt.show()
