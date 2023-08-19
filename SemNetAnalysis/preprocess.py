"""
__________________________________________________________________________________________________________________________________________
File: preprocess.py
Author: Sahana Mahesh
Email: sahana.sm61@gmail.com
Date: July 8th, 2023
Last Modified: August 1st, 2023

Description: This script is part of a COGS 402 project (UBC 2023S tearm 1 & 2) conducted at the Cognitive Neuroscience of Thought Laboratory 
             under the supervision of Dr. Kalina Christoff (PI) and Jen Burrell. It takes a .txt file, calculates the texts word count, 
             removes punctuations, numbers, dates, stopwords, performs lemmatizations and creates individual tokens (words). 
             It returns a list of words.
__________________________________________________________________________________________________________________________________________
"""
import nltk
from nltk.stem import PorterStemmer
from nltk.tokenize import word_tokenize, MWETokenizer
from nltk.corpus import stopwords
from collections import Counter
import string
import re

class Preprocessor:

    def __init__(self, file, ngrams=[]):
        self.file = file
        self.ngrams = ngrams
        self.text = self.read_text()
        self.words = []
        self.word_count = 0

    # read_text(file) reades a given .txt file
    def read_text(self):
        with open(self.file, 'r') as f:
            text = f.read()
        return text
    
    # Count the number of words before cleaning
    def get_word_count(self):
        lines = self.text.split()
        self.word_count += len(lines)
        return self.word_count


    def lowercase_text(self):
        self.text = self.text.lower()
        return self.text

    # remove_punctuation(self) function takes a text file, removes punctuations, 
    # and returns a text file without punctuations
    def remove_punctuation(self):
        self.text = self.text.translate(str.maketrans('', '', string.punctuation))
        return self.text
    
    # remove_numbers_and_dates(self) function takes a text file, removes numbers and dates, 
    # and returns a text file
    def remove_numbers_and_dates(self):
        self.text = re.sub(r'\b\d+\b', '', self.text)
        self.text = re.sub(r'\b(\d{1,2}[-/]\d{1,2}[-/]\d{2,4})\b', '', self.text)
        return self.text

    # tokenize_word(self) function segments the text into individual word units 
    # and returns a list of words
    def tokenize_word(self):
        tokenizer = MWETokenizer(self.ngrams)
        self.words = tokenizer.tokenize(word_tokenize(self.text))
        self.words = [word.replace('_', ' ') for word in self.words]
        return self.words

    # remove_stopwords(self) function takes a list of words, removes stopwords from the list,
    # and returns a list of words without stopwords
    def remove_stopwords(self):
        stop_words = set(stopwords.words('english'))
        filler_words = ['like', 'um', 'uh', 'er', 'ah', 'yeah', 'na', 'im', 'oh', 'hmm']
        stop_words.update(filler_words)
        self.words = [word for word in self.words if word not in stop_words]
        return self.words

    # stem_and_lemmatiza(self) function takes a list of words and returns a list of word 
    # (all words in words are reduced to root form), performs stemming to remove suffixes and 
    # lemmatization to remove synonyms
    def stem_and_lemmatize(self):
        #porter = PorterStemmer()
        lemmatizer = nltk.WordNetLemmatizer()
        #self.words = [porter.stem(word) for word in self.words] # did we still want to include this?
        self.words = [lemmatizer.lemmatize(word) for word in self.words] 
        return self.words
         
    def preprocess(self):
        # Performs all preprocessing steps in sequence
        self.lowercase_text()
        self.remove_punctuation()
        self.remove_numbers_and_dates()
        self.tokenize_word()
        self.remove_stopwords()
        self.stem_and_lemmatize()
        return self.words
        
    def write_to_file(self, filename):
        #preprocessed_text = ' '.join(self.preprocess())
        with open(filename, 'w') as f:
            f.write(' '.join(self.words))
        #print(f"Preprocessed data written to {filename}")