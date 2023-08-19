"""
__________________________________________________________________________________________________________________________________________
File: cooccurrence_network.py
Author: Sahana Mahesh
Email: sahana.sm61@gmail.com
Date: July 10th, 2023
Last Modified: August 1st, 2023

Description: This script is part of a COGS 402 project (UBC 2023S tearm 1 & 2) conducted at the Cognitive Neuroscience of Thought Laboratory 
             under the supervision of Dr. Kalina Christoff (PI) and Jen Burrell. This script estimates a cooccurrence semantic network 
             with a window size of 10 (window size = 5 by default). The script take a list of tokens, each token is represented as a node
             in the network. The edge weights in the network are estimated based on the Inverse Distance Weighting (IDW) scheme.
__________________________________________________________________________________________________________________________________________
"""
import numpy as np
import networkx as nx
from collections import Counter
from scipy.sparse import csr_matrix


class CooccurrenceNetworkBuilder:

    # Initialize CooccurrenceNetworkBuilder class
    def __init__(self, tokens, window_size=10):
        
        self.tokens = tokens
        self.window_size = window_size

    def build_cooccurrence_network(self):

        # Count frequency of each word
        counter = Counter(self.tokens)

        # Map each word to a unique id
        vocab = {word: i for i, word in enumerate(counter.keys())}

    # Get the total number of unique words
        vocab_size = len(vocab)

    # Prepare an empty matrix to store cooccurrences
        cooccurrence_matrix = np.zeros((vocab_size, vocab_size), dtype=np.float32)

    # Iterate over each word in the text
        for i, token in enumerate(self.tokens):
        # Check each word within the window around the current word
            for offset in range(1, self.window_size+1):
            # Check the word at the position offset to the left
                if i - offset >= 0:
                # Add weight based on distance
                    cooccurrence_matrix[vocab[token]][vocab[self.tokens[i - offset]]] += 1 / offset
            # Check the word at the position offset to the right
                if i + offset < len(self.tokens):
                # Add weight based on distance
                    cooccurrence_matrix[vocab[token]][vocab[self.tokens[i + offset]]] += 1 / offset

    # Convert cooccurrence matrix to a sparse format (efficiency)
        cooccurrence_matrix = csr_matrix(cooccurrence_matrix)

    # Convert CSR matrix to a NetworkX graph
        network = nx.from_numpy_array(cooccurrence_matrix.toarray())
        

    # Create a mapping from id to word
        id_to_word = {i: word for word, i in vocab.items()}

    # Relabel nodes with their corresponding words
        network = nx.relabel_nodes(network, id_to_word)
        #print(network.edges(data=True))
        return network
