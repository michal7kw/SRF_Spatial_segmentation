import pandas as pd

# Read the CSV file
df = pd.read_csv('detected_transcripts.csv')

# Set all values in the 'global_z' column to 0.0
df['global_z'] = 0.0

# Save the modified dataframe back to the CSV file
df.to_csv('detected_transcripts.csv', index=False)