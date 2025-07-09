import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file into a DataFrame
df = pd.read_csv('rmsd-over-time-backbone.csv')

# Strip any leading/trailing whitespace from column names
df.columns = df.columns.str.strip()

# Rename 'Time (ps)' to 'Time' for simplicity
df.rename(columns={'Time (ps)': 'Time'}, inplace=True)

# Define the list of replicate columns
rmsd_cols = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10']

# Melt the DataFrame to long format for easier plotting with Seaborn
melted_df = df.melt(id_vars=['Time'], value_vars=rmsd_cols, 
                    var_name='Replicate', value_name='RMSD')

# Plot RMSD over time for all replicates
plt.figure(figsize=(10, 6))
sns.lineplot(data=melted_df, x='Time', y='RMSD', hue='Replicate')
plt.title('Backbone RMSD for Each Replicate')
plt.xlabel('Time (ps)')
plt.ylabel('RMSD (nm)')
plt.show()

# Compute mean and standard deviation across replicates at each time point
df['Mean'] = df[rmsd_cols].mean(axis=1)
df['Std'] = df[rmsd_cols].std(axis=1)

# Plot mean RMSD with standard deviation as a shaded band
plt.figure(figsize=(10, 6))
plt.plot(df['Time'], df['Mean'], label='Mean', color='black')
plt.fill_between(df['Time'], df['Mean'] - df['Std'], df['Mean'] + df['Std'], 
                 alpha=0.2, color='gray')
plt.legend()
plt.title('Backbone RMSD: Mean with Standard Deviation')
plt.xlabel('Time (ps)')
plt.ylabel('RMSD (nm)')
plt.show()

# Define the equilibrated phase (last 50% of the trajectory)
total_time = df['Time'].max()
equil_start = total_time / 2  # Starts at 75 ps for 150 ps total

# Filter the DataFrame for the equilibrated phase
equil_df = df[df['Time'] > equil_start]

# Compute mean and standard deviation for each replicate in the equilibrated phase
equil_means = equil_df[rmsd_cols].mean()
equil_stds = equil_df[rmsd_cols].std()

# Print equilibrated phase statistics
print("Equilibrated Backbone RMSD:")
for rep in rmsd_cols:
    print(f"{rep}: Mean = {equil_means[rep]:.3f} nm, Std = {equil_stds[rep]:.3f} nm")

# Plot bar plot of mean equilibrated RMSD with error bars
plt.figure(figsize=(10, 6))
equil_means.plot(kind='bar', yerr=equil_stds, capsize=5)
plt.title('Mean Equilibrated Backbone RMSD Across Replicates')
plt.ylabel('RMSD (nm)')
plt.show()

# Compute correlation matrix of RMSD across replicates
corr_matrix = df[rmsd_cols].corr()

# Plot heatmap of the correlation matrix
plt.figure(figsize=(8, 6))
sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
plt.title('Correlation Matrix of Backbone RMSD Across Replicates')
plt.show()
