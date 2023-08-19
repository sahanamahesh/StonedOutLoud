"""
__________________________________________________________________________________________________________________________________________
File: preprocess.py
Author: Sahana Mahesh
Email: sahana.sm61@gmail.com
Date: July 10th, 2023
Last Modified: August 1st, 2023

Description: This script is part of a COGS 402 project (UBC 2023S tearm 1 & 2) conducted at the Cognitive Neuroscience of Thought Laboratory 
             under the supervision of Dr. Kalina Christoff (PI) and Jen Burrell. This script performs various tests to assure that the 
             methods in the Preprocessor class are functioning as required.
__________________________________________________________________________________________________________________________________________
"""
import unittest
from nltk.corpus import stopwords
from nltk.stem import PorterStemmer, WordNetLemmatizer
from preprocess import Preprocessor 


class TestPreprocessor(unittest.TestCase):

    def setUp(self):
        self.preprocessor = Preprocessor('TestData copy/sub-test-002_ses-1_task-TAP.txt')

    def test_read_text(self):
        # read the text from file
        self.preprocessor.read_text()
        self.assertIsNotNone(self.preprocessor.text)

    def test_lowercase_text(self):
        self.preprocessor.text = "TEST"
        self.preprocessor.lowercase_text()
        self.assertEqual(self.preprocessor.text, "test")

    def test_remove_punctuation(self):
        self.preprocessor.text = "hello, world!"
        self.preprocessor.remove_punctuation()
        self.assertEqual(self.preprocessor.text, "hello world")

    def test_remove_numbers_and_dates(self):
        self.preprocessor.text = "text with 123 and date 20-12-2020"
        self.preprocessor.remove_numbers_and_dates()
        self.preprocessor.remove_punctuation()
        self.assertEqual(self.preprocessor.text, "text with  and date ")

    def test_tokenize_word(self):
        self.preprocessor.text = "This is a sentence"
        self.preprocessor.tokenize_word()
        self.assertEqual(self.preprocessor.words, ['This', 'is', 'a', 'sentence'])

    def test_remove_stopwords(self):
        self.preprocessor.words = ['This', 'is', 'a', 'sentence']
        self.preprocessor.remove_stopwords()
        self.assertEqual(self.preprocessor.words, ['This', 'sentence'])

    def test_stem_and_lemmatize(self):
        self.preprocessor.words = ['played', 'running']
        self.preprocessor.stem_and_lemmatize()
        self.assertEqual(self.preprocessor.words, ['play', 'run'])

    def test_preprocess(self):
        words = self.preprocessor.preprocess()
        self.assertIsNotNone(words)

if __name__ == '__main__':
    unittest.main()
