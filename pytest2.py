import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file into a DataFrame with the correct delimiter
df = pd.read_csv('rg-over-time-protein.csv', sep='\s+')

# Verify column names
print(df.columns)

# Define the list of replicate columns
rg_cols = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10']

# Melt the DataFrame to long format for easier plotting with Seaborn
melted_df = df.melt(id_vars=['Time (ps)'], value_vars=rg_cols, 
                    var_name='Replicate', value_name='Rg')

# Plot Rg over time for all replicates
plt.figure(figsize=(10, 6))
sns.lineplot(data=melted_df, x='Time (ps)', y='Rg', hue='Replicate')
plt.title('Radius of Gyration (Rg) for Each Replicate')
plt.xlabel('Time (ps)')
plt.ylabel('Rg (nm)')
plt.show()

# Compute mean and standard deviation across replicates at each time point
df['Mean'] = df[rg_cols].mean(axis=1)
df['Std'] = df[rg_cols].std(axis=1)

# Plot mean Rg with standard deviation as a shaded band
plt.figure(figsize=(10, 6))
plt.plot(df['Time (ps)'], df['Mean'], label='Mean', color='black')
plt.fill_between(df['Time (ps)'], df['Mean'] - df['Std'], df['Mean'] + df['Std'], 
                 alpha=0.2, color='gray')
plt.legend()
plt.title('Radius of Gyration (Rg): Mean with Standard Deviation')
plt.xlabel('Time (ps)')
plt.ylabel('Rg (nm)')
plt.show()

# Define the equilibrated phase (last 50% of the trajectory)
total_time = df['Time (ps)'].max()
equil_start = total_time / 2  # Starts at 75 ps for 150 ps total

# Filter the DataFrame for the equilibrated phase
equil_df = df[df['Time (ps)'] > equil_start]

# Compute mean and standard deviation for each replicate in the equilibrated phase
equil_means = equil_df[rg_cols].mean()
equil_stds = equil_df[rg_cols].std()

# Print equilibrated phase statistics
print("Equilibrated Radius of Gyration (Rg):")
for rep in rg_cols:
    print(f"{rep}: Mean = {equil_means[rep]:.3f} nm, Std = {equil_stds[rep]:.3f} nm")

# Plot bar plot of mean equilibrated Rg with error bars
plt.figure(figsize=(10, 6))
equil_means.plot(kind='bar', yerr=equil_stds, capsize=5)
plt.title('Mean Equilibrated Radius of Gyration (Rg) Across Replicates')
plt.ylabel('Rg (nm)')
plt.show()

# Compute correlation matrix of Rg across replicates
corr_matrix = df[rg_cols].corr()

# Plot heatmap of the correlation matrix
plt.figure(figsize=(8, 6))
sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
plt.title('Correlation Matrix of Radius of Gyration (Rg) Across Replicates')
plt.show()
