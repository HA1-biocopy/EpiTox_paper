"""
Clean Publication-Quality Radar Plot Generator
Creates the radar plot with informative value labels exactly as shown in the PDF
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from math import pi

# Read data
dir = "../../data/"
df = pd.read_excel('../../data/experimetnal_data_table.xlsx')

# Calculate statistics for each unique peptide ID
id_stats = []
for peptide_id in df['id'].unique():
    id_data = df[df['id'] == peptide_id]
    stats = {
        'id': peptide_id,
        'n_observations': len(id_data),
        'n_diseases': id_data['disease'].dropna().nunique(),
        'n_hla_alleles': id_data['hla_allele'].nunique(),
        'n_hla_classes': id_data['hla_class'].nunique(),
        'n_studies': id_data['study_ref'].nunique(),
        'n_methods': id_data['Experimental_method'].nunique(),
        'n_binding_types': id_data['Binding'].nunique()
    }
    id_stats.append(stats)

id_df = pd.DataFrame(id_stats)

# Calculate average across all IDs
avg_stats = {
    'Total\nObservations': id_df['n_observations'].mean(),
    'Diseases': id_df['n_diseases'].mean(),
    'HLA\nAlleles': id_df['n_hla_alleles'].mean(),
    'HLA\nClasses': id_df['n_hla_classes'].mean(),
    'Studies': id_df['n_studies'].mean(),
    'Experimental\nMethods': id_df['n_methods'].mean(),
    'Binding\nTypes': id_df['n_binding_types'].mean()
}

# Print statistics
print("="*60)
print("AVERAGE DATA COVERAGE PER UNIQUE PEPTIDE ID")
print("="*60)
for key, val in avg_stats.items():
    key_clean = key.replace('\n', ' ')
    print(f"{key_clean:.<30} {val:.2f}")
print(f"\nTotal unique peptides: {len(id_df)}")
print(f"Total observations: {len(df)}")
print("="*60)

# Create the clean publication-quality radar plot
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))

# Prepare data
categories = list(avg_stats.keys())
values = list(avg_stats.values())

# Number of variables
N = len(categories)

# Compute angle for each axis
angles = [n / float(N) * 2 * pi for n in range(N)]
values += values[:1]  # Complete the circle
angles += angles[:1]

# Plot data
ax.plot(angles, values, 'o-', linewidth=3, color='#2E86AB', markersize=10)
ax.fill(angles, values, alpha=0.3, color='#2E86AB')

# Fix axis to go in the right order and start at 12 o'clock
ax.set_theta_offset(pi / 2)
ax.set_theta_direction(-1)

# Draw axis lines for each angle and label
ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories, size=12, weight='bold')

# Set y-axis limits with some padding for labels
max_value = max(values)
ax.set_ylim(0, max_value * 1.3)

# Add gridlines with custom values
gridline_values = [2, 4, 6, 8]
ax.set_yticks(gridline_values)
ax.set_yticklabels([str(x) for x in gridline_values], size=10)

# Grid styling
ax.grid(True, linewidth=1.5, alpha=0.3)

# Add value labels on the plot with nice formatting
for angle, value, label in zip(angles[:-1], values[:-1], categories):
    # Position labels slightly outside the data points
    ax.text(angle, value * 1.15, f'{value:.2f}', 
            ha='center', va='center', size=11, weight='bold',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', 
                     edgecolor='#2E86AB', linewidth=2, alpha=0.9))

# Title with subtitle
ax.set_title('Average data dimensions per Uniprot-Peptide pair\n' + 
             '258 unique peptides | 2198 total observations\n' ,
            size=14, weight='bold', pad=35, loc='center')
# Adjust layout
plt.tight_layout()

# Save in multiple formats
#plt.savefig('radar_plot_clean_publication.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(dir+'radar_plot_clean.pdf', dpi=300, bbox_inches='tight', facecolor='white')

print("\n✓ Files saved:")
print("  - radar_plot_clean_publication.png")
print("  - radar_plot_clean_publication.pdf")

# Also create a version with additional context
fig2, ax2 = plt.subplots(figsize=(12, 10), subplot_kw=dict(projection='polar'))

# Plot with same data
ax2.plot(angles, values, 'o-', linewidth=3.5, color='#2E86AB', markersize=12)
ax2.fill(angles, values, alpha=0.25, color='#2E86AB')

ax2.set_theta_offset(pi / 2)
ax2.set_theta_direction(-1)
ax2.set_xticks(angles[:-1])
ax2.set_xticklabels(categories, size=13, weight='bold')
ax2.set_ylim(0, max_value * 1.4)
ax2.set_yticks(gridline_values)
ax2.set_yticklabels([str(x) for x in gridline_values], size=11)
ax2.grid(True, linewidth=1.5, alpha=0.3)

# Add value labels with percentage of total where applicable
for i, (angle, value, label) in enumerate(zip(angles[:-1], values[:-1], categories)):
    # Calculate percentage of total for relevant dimensions
    if 'Diseases' in label:
        pct = (value / df['disease'].dropna().nunique()) * 100
        label_text = f'{value:.2f}\n({pct:.1f}% of total)'
    elif 'HLA\nAlleles' in label:
        pct = (value / df['hla_allele'].nunique()) * 100
        label_text = f'{value:.2f}\n({pct:.1f}% of total)'
    elif 'Studies' in label:
        pct = (value / df['study_ref'].nunique()) * 100
        label_text = f'{value:.2f}\n({pct:.1f}% of total)'
    else:
        label_text = f'{value:.2f}'
    
    ax2.text(angle, value * 1.18, label_text, 
            ha='center', va='center', size=10, weight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', 
                     edgecolor='#2E86AB', linewidth=2, alpha=0.95))

ax2.set_title('Average data dimensions per Uniprot-Peptide pair\n' + 
             '258 unique peptides | 2198 total observations\n' ,
            size=14, weight='bold', pad=35, loc='center')

plt.tight_layout()
#plt.savefig('radar_plot_with_percentages.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(dir+'radar_plot_with_percentages.pdf', dpi=300, bbox_inches='tight', facecolor='white')

#print("  - radar_plot_with_percentages.png")
print("  - radar_plot_with_percentages.pdf")

# Create version with comparison (min, avg, max)
fig3, ax3 = plt.subplots(figsize=(12, 10), subplot_kw=dict(projection='polar'))

# Calculate min and max
min_stats = {
    'Total\nObservations': id_df['n_observations'].min(),
    'Diseases': id_df['n_diseases'].min(),
    'HLA\nAlleles': id_df['n_hla_alleles'].min(),
    'HLA\nClasses': id_df['n_hla_classes'].min(),
    'Studies': id_df['n_studies'].min(),
    'Experimental\nMethods': id_df['n_methods'].min(),
    'Binding\nTypes': id_df['n_binding_types'].min()
}

max_stats = {
    'Total\nObservations': id_df['n_observations'].max(),
    'Diseases': id_df['n_diseases'].max(),
    'HLA\nAlleles': id_df['n_hla_alleles'].max(),
    'HLA\nClasses': id_df['n_hla_classes'].max(),
    'Studies': id_df['n_studies'].max(),
    'Experimental\nMethods': id_df['n_methods'].max(),
    'Binding\nTypes': id_df['n_binding_types'].max()
}

min_values = list(min_stats.values()) + [list(min_stats.values())[0]]
max_values = list(max_stats.values()) + [list(max_stats.values())[0]]

# Plot all three
ax3.plot(angles, min_values, 'o-', linewidth=2, color='#A23B72', 
        label='Minimum', alpha=0.7, markersize=8)
ax3.fill(angles, min_values, alpha=0.15, color='#A23B72')

ax3.plot(angles, values, 'o-', linewidth=3.5, color='#2E86AB', 
        label='Average', alpha=0.9, markersize=12)
ax3.fill(angles, values, alpha=0.25, color='#2E86AB')

ax3.plot(angles, max_values, 'o-', linewidth=2, color='#F18F01', 
        label='Maximum', alpha=0.7, markersize=8)
ax3.fill(angles, max_values, alpha=0.15, color='#F18F01')

ax3.set_theta_offset(pi / 2)
ax3.set_theta_direction(-1)
ax3.set_xticks(angles[:-1])
ax3.set_xticklabels(categories, size=13, weight='bold')
ax3.set_ylim(0, max(max_values) * 1.1)
ax3.grid(True, linewidth=1.5, alpha=0.3)

# Add legend
ax3.legend(loc='upper right', bbox_to_anchor=(1.25, 1.1), fontsize=12, frameon=True, 
          fancybox=True, shadow=True)

ax3.set_title('Data Dimension Range Across All Peptide IDs\n(Minimum, Average, Maximum)',
            size=14, weight='bold', pad=35)

plt.tight_layout()
#plt.savefig('radar_plot_comparison.png', dpi=300, bbox_inches='tight', facecolor='white')
#plt.savefig('radar_plot_comparison.pdf', dpi=300, bbox_inches='tight', facecolor='white')

# print("  - radar_plot_comparison.png")
# print("  - radar_plot_comparison.pdf")

print("\n" + "="*60)
print("ALL RADAR PLOTS GENERATED SUCCESSFULLY!")
print("="*60)
