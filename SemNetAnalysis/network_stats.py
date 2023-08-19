"""
__________________________________________________________________________________________________________________________________________
File: network_stats.py
Author: Sahana Mahesh
Co-Author: Jen Burrell
Email: sahana.sm61@gmail.com
Date: July 20th, 2023
Last Modified: August 1st, 2023

Description: This script is part of a COGS 402 project (UBC 2023S tearm 1 & 2) conducted at the Cognitive Neuroscience of Thought Laboratory 
             under the supervision of Dr. Kalina Christoff (PI) and Jen Burrell. It takes a graph with nodes and edges to calculate and 
             returns a data frame containing the following network measures: number of nodes, Number of edges, average degree, average 
             betweenness centrality, average clustering coefficient, average path length, modularity, number of communities, and word count.
__________________________________________________________________________________________________________________________________________
"""
import networkx as nx
import pandas as pd
from statistics import mean
from community import community_louvain


class NetworkStats:
    def __init__(self, network):
        # Initalize variables
        self.network = network
        self.stats_df = pd.DataFrame()

    def calculate_stats(self):

        stats = {}

        G = nx.Graph(self.network)

        number_of_nodes = len(G.nodes())
        stats['Nodes'] = number_of_nodes

        number_of_edges = len(G.edges())
        stats['Edges'] = number_of_edges
        
        degree_centrality = dict(nx.degree(self.network, weight = 'weight'))
        stats['degree_centrality'] = mean(degree_centrality.values())
        
        betweenness_centrality = nx.betweenness_centrality(self.network, weight='weight')
        stats['betweenness_centrality'] = mean(betweenness_centrality.values())

        # closeness_centrality = nx.closeness_centrality(self.network, distance='weight')
        # stats['closeness_centrality'] = mean(closeness_centrality.values())

        clustering_coefficient = nx.clustering(self.network, nodes=None, weight="weight")
        stats['clustering_coefficient'] = mean(clustering_coefficient.values())

        stats['average_path_length'] = nx.average_shortest_path_length(self.network, weight='weight')

        # eigenvector_centrality = nx.eigenvector_centrality_numpy(self.network, weight='weight')
        # stats['eigenvector_centrality'] = mean(eigenvector_centrality.values())

        partition = community_louvain.best_partition(self.network, resolution=1)
        # modularity = community_louvain.modularity(partition, self.network)
        stats['modularity'] = community_louvain.modularity(community_louvain.best_partition(self.network, resolution=1), self.network)

        communities = len(set(partition.values()))

        stats['communities'] = communities
        
        self.stats_df = pd.DataFrame([stats])

        return self.stats_df

    
    # def calculate_individual_stats(self):
    #     stats = pd.DataFrame()

    #     degree = nx.degree(self.network, weight='weight')
    #     stats['degree'] = degree

    #     degree_centrality = dict(nx.degree(self.network, weight = 'weight'))
    #     stats['degree_centrality'] = degree_centrality
        
    #     betweenness_centrality = nx.betweenness_centrality(self.network, weight='weight')
    #     stats['betweenness_centrality'] = betweenness_centrality

    #     closeness_centrality = nx.closeness_centrality(self.network, distance='weight')
    #     stats['closeness_centrality'] = closeness_centrality

    #     clustering_coefficient = nx.clustering(self.network, nodes=None, weight="weight")
    #     stats['clustering_coefficient'] = clustering_coefficient

    #     eigenvector_centrality = nx.eigenvector_centrality_numpy(self.network, weight='weight')
    #     stats['eigenvector_centrality'] = eigenvector_centrality

    #     return stats