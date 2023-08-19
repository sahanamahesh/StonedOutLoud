"""
__________________________________________________________________________________________________________________________________________
File: generate_plots.py
Author: Sahana Mahesh
Email: sahana.sm61@gmail.com
Date: August 8th, 2023
Last Modified: August 10th, 2023

Description: This script is part of a COGS 402 project (UBC 2023S tearm 1 & 2) conducted at the Cognitive Neuroscience of Thought Laboratory 
             under the supervision of Dr. Kalina Christoff (PI) and Jen Burrell. It takes an excel file and generates box plots with network 
             measures on the x axis: number of nodes, Number of edges, average degree, average betweenness centrality, average clustering  
             coefficient, average path length, modularity, number of communities, and word count. The y axis has normalized values [1, 0].
__________________________________________________________________________________________________________________________________________
"""
import pandas as pd
import matplotlib.pyplot as plt

class GeneratePlots:
    def __init__(self, filepath):
        self.data = pd.read_excel(filepath)
        self.features = ['Word Count', 'Nodes', 'Edges', 'Degree Centrality', 'Betweenness Centrality', 
                'Clustering Coefficient', 'Average Path Length',
                'Modularity', 'Communities']
        self._normalize_features()

        # Define groups
        self.groups = {
            'Smoking Vs Sober': {
                'Smoking': self.data[self.data['Smoking'] == 1],
                'Sober': self.data[self.data['Smoking'] == 0]
            },
            'Session 1 Smoking Vs Session 2 Sober': {
                'Ses-1 Smoking': self.data[(self.data['Smoking'] == 1) & (self.data['Session'] == 1)],
                'Ses-2 Sober': self.data[(self.data['Smoking'] == 0) & (self.data['Session'] == 2)]
            },
            'Session 2 Smoking Vs Session 1 Sober': {
                'Ses-2 Smoking': self.data[(self.data['Smoking'] == 1) & (self.data['Session'] == 2)],
                'Ses-1 Sober': self.data[(self.data['Smoking'] == 0) & (self.data['Session'] == 1)]
            },
            'Session 1 Smoking Vs Session 2 Smoking': {
                'Ses-1 Smoking': self.data[(self.data['Smoking'] == 1) & (self.data['Session'] == 1)],
                'Ses-2 Smoking': self.data[(self.data['Smoking'] == 1) & (self.data['Session'] == 2)]
            },
            'Session 1 Sober Vs Session 2 Sober': {
                'Ses-1 Smoking': self.data[(self.data['Smoking'] == 0) & (self.data['Session'] == 1)],
                'Ses-2 Smoking': self.data[(self.data['Smoking'] == 0) & (self.data['Session'] == 2)]
            }
        }

    def _normalize_features(self):
        # Normalize the feature values
        for feature in self.features:
            self.data[feature] = (self.data[feature] - self.data[feature].min()) / (self.data[feature].max() - self.data[feature].min())

    def plot_boxplot(self, group_key, save_path=None):
        group_data = self.groups[group_key]

        plot_data = {label: subset[self.features] for label, subset in group_data.items()}
        
        fig, ax = plt.subplots(figsize=(35, 10))
        
        # Adjust spacing
        inner_spacing = 1
        outer_spacing = 1  

        # Iterate over each feature to create side-by-side box plots
        for i, feature in enumerate(self.features):
            data_to_plot = [group[feature].values for group in plot_data.values()]
            positions = range(i * (len(plot_data) * inner_spacing + outer_spacing), 
                              i * (len(plot_data) * inner_spacing + outer_spacing) + len(plot_data) * inner_spacing, 
                              inner_spacing)
            box = ax.boxplot(data_to_plot, positions=positions, widths=0.8, patch_artist=True, medianprops=dict(color='black', linewidth=2))
            
            for patch, color, hatch in zip(box['boxes'], ['white', 'white'], ['/', '.']):
                patch.set_facecolor(color)
                patch.set_hatch(hatch)
                patch.set_edgecolor('black')
                
        ax.set_xticks([(i * (len(plot_data) * inner_spacing + outer_spacing) + (len(plot_data)-1) * 0.5 * inner_spacing) for i in range(len(self.features))])
        ax.set_xticklabels(self.features, fontsize=14)
        ax.set_title('Distribution of Network Measures for ' + group_key, fontsize=20)
        ax.set_xlabel('Network Measures', fontsize=16)
        ax.set_ylabel('Normalized Value', fontsize=16)
        
        # Adding a custom legend
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor='white', edgecolor='black', hatch=hatch, label=label) for label, hatch in zip(plot_data.keys(), ['/', '.'])]
        ax.legend(handles=legend_elements, title='Group', loc='center left', bbox_to_anchor=(1, 0.5), fontsize=14)
        
        # Adjust x-axis limits
        ax.set_xlim(-2, (len(self.features) - 1) * (len(plot_data) * inner_spacing + outer_spacing) + len(plot_data) * inner_spacing + 2)
        
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path + group_key + '.png')

        #plt.show()

if __name__ == "__main__":
    visualizer = GeneratePlots("Derivatives/all_subs/SemNetAnalysis/Plots/network_stats_plot.xlsx")

    # Plot and save for each group
    groups_to_plot = ['Smoking Vs Sober', 'Session 1 Smoking Vs Session 2 Sober', 'Session 2 Smoking Vs Session 1 Sober', 'Session 1 Smoking Vs Session 2 Smoking', 'Session 1 Sober Vs Session 2 Sober']
    save_path = 'Derivatives/all_subs/SemNetAnalysis/Plots/'

    for group in groups_to_plot:
        visualizer.plot_boxplot(group, save_path)

