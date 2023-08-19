"""
__________________________________________________________________________________________________________________________________________
File: analysis.py
Author: Sahana Mahesh
Email: sahana.sm61@gmail.com
Date: July 21st, 2023
Last Modified: August 10th, 2023

Description: This script is part of a COGS 402 project (UBC 2023S tearm 1 & 2) conducted at the Cognitive Neuroscience of Thought Laboratory 
             under the supervision of Dr. Kalina Christoff (PI) and Jen Burrell. It takes a data frame and a list of features containing 
             label name. It calculates descriptive statistis, performs the Wilcoxon signed-rank test, and the Mann-Whitney U test, corrects 
             for multiple comparisons using the Bonferroni-Holm method and saves the results in individual excel files.
__________________________________________________________________________________________________________________________________________
"""
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon, mannwhitneyu, norm
from statsmodels.sandbox.stats.multicomp import multipletests
import os
import math

class Analysis:
    # Class initialization
    def __init__(self, data, features):
        self.data = data
        self.features = features
        self.alpha = 0.05

    # Function to run the Wilcoxon test on two groups
    def run_wilcoxon_test(self, group1, group2):
        result = wilcoxon(group1, group2)
        # N = math.sqrt(len(group1))
        # effect_size = result.statistic / N 
        n = len(group1)
        rank_biserial_r = 1 - (2 * result.statistic) / (n * (n + 1))
        effect_size = rank_biserial_r

        return result.statistic, result.pvalue, effect_size

    # Function to run the Mann Whitney U test on two groups
    def run_mann_whitney_test(self, group1, group2):
        result = mannwhitneyu(group1, group2)
        # Z = abs(norm.ppf(result.pvalue / 2))
        # N = len(group1) + len(group2)
        # effect_size = Z / np.sqrt(N)

        n1 = len(group1)
        n2 = len(group2)
        mean_U = n1 * n2 / 2
        std_U = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
        
        Z = (result.statistic - mean_U) / std_U
        N = n1 + n2
        effect_size = Z / np.sqrt(N)

        return result.statistic, result.pvalue, Z, effect_size
    
    # Function to calculate mean, median, sd, and iqr for each group
    def calculate_descriptive_stats(self, group):

        mean = group.mean()
        median = group.median()
        sd = group.std()
        iqr = np.percentile(group, 75) - np.percentile(group, 25)

        return(mean, median, sd, iqr)

    # Function to analyze the data for smoke and sober groups
    def run_smoke_sober_analysis(self):
        results_smoke_sober = {'Feature': [], 'Statistic': [], 'p-value': [], 'Effect size': [], 
                               'mean_smoke':[], 'mean_sober':[], 'median_smoke':[], 'median_sober':[],
                                'sd_smoke': [], 'sd_sober': [], 'iqr_smoke': [], 'iqr_sober': []}

        smoke = self.data[(self.data['Smoking'] == 1)]
        sober = self.data[(self.data['Smoking'] == 0)]

        for feature in self.features:
            statistic, p, r = self.run_wilcoxon_test(smoke[feature], sober[feature])
            mean_smoke, median_smoke, sd_smoke, iqr_smoke = self.calculate_descriptive_stats(smoke[feature])
            mean_sober, median_sober, sd_sober, iqr_sober = self.calculate_descriptive_stats(sober[feature])

            results_smoke_sober['Feature'].append(feature)
            results_smoke_sober['Statistic'].append(statistic)
            results_smoke_sober['p-value'].append(p)
            results_smoke_sober['Effect size'].append(r)
            results_smoke_sober['mean_smoke'].append(mean_smoke)
            results_smoke_sober['median_smoke'].append(median_smoke)
            results_smoke_sober['sd_smoke'].append(sd_smoke)
            results_smoke_sober['iqr_smoke'].append(iqr_smoke)
            results_smoke_sober['mean_sober'].append(mean_sober)
            results_smoke_sober['median_sober'].append(median_sober)
            results_smoke_sober['sd_sober'].append(sd_sober)
            results_smoke_sober['iqr_sober'].append(iqr_sober)

        reject, pvals_corrected, _, _ = multipletests(results_smoke_sober['p-value'], method='holm', alpha=self.alpha)
        results_smoke_sober['p-value_corrected'] = pvals_corrected
        results_smoke_sober['reject'] = reject

        results_smoke_sober = pd.DataFrame(results_smoke_sober)
        results_smoke_sober.to_excel(os.path.join('Derivatives/all_subs/SemNetAnalysis/group_analyses_results', 'sober_vs_cannabis_group_level.xlsx'), index=False)

    # Function to analyze the data for TAP and FAST groups
    def run_TAP_FAST_analysis(self):
        results_tap_fast = {'Feature': [], 'tap_Statistic': [], 'tap_p-value': [], 'tap Effect size': [], 
                            'fast_Statistic': [], 'fast_p-value': [], 'fast Effect size': [], 
                            'smoke_Statistic': [], 'smoke_p-value': [], 'smoke_Z': [], 'smoke Effect size': [], 
                            'sober_Statistic': [], 'sober_p-value': [], 'sober_Z': [], 'sober Effect size': []}
        
        results_tap_fast_d = {'Feature': [], 'tap_mean_smoke': [], 'tap_mean_sober': [], 'fast_mean_smoke': [], 'fast_mean_sober': [],
                              'tap_median_smoke': [], 'tap_median_sober': [], 'fast_median_smoke': [], 'fast_median_sober': [],
                              'tap_sd_smoke': [], 'tap_sd_sober': [], 'fast_sd_smoke': [], 'fast_sd_sober': [],
                              'tap_iqr_smoke': [], 'tap_iqr_sober': [], 'fast_iqr_smoke': [], 'fast_iqr_sober': []}

        tap_smoke = self.data[(self.data['Smoking'] == 1) & (self.data['TAP'] == 1)]
        tap_sober = self.data[(self.data['Smoking'] == 0) & (self.data['TAP'] == 1)]
        fast_smoke = self.data[(self.data['Smoking'] == 1) & (self.data['TAP'] == 0)]
        fast_sober = self.data[(self.data['Smoking'] == 0) & (self.data['TAP'] == 0)]

        for feature in self.features:
            mean_smoke, median_smoke, sd_smoke, iqr_smoke = self.calculate_descriptive_stats(tap_smoke[feature])
            mean_sober, median_sober, sd_sober, iqr_sober = self.calculate_descriptive_stats(tap_sober[feature])

            results_tap_fast_d['Feature'].append(feature)
            results_tap_fast_d['tap_mean_smoke'].append(mean_smoke)
            results_tap_fast_d['tap_mean_sober'].append(mean_sober)
            results_tap_fast_d['tap_median_smoke'].append(median_smoke)
            results_tap_fast_d['tap_median_sober'].append(median_sober)
            results_tap_fast_d['tap_sd_smoke'].append(sd_smoke)
            results_tap_fast_d['tap_sd_sober'].append(sd_sober)
            results_tap_fast_d['tap_iqr_smoke'].append(iqr_smoke)
            results_tap_fast_d['tap_iqr_sober'].append(iqr_sober)

            mean_smoke, median_smoke, sd_smoke, iqr_smoke = self.calculate_descriptive_stats(fast_smoke[feature])
            mean_sober, median_sober, sd_sober, iqr_sober = self.calculate_descriptive_stats(fast_sober[feature])

            results_tap_fast_d['fast_mean_smoke'].append(mean_smoke)
            results_tap_fast_d['fast_mean_sober'].append(mean_sober)
            results_tap_fast_d['fast_median_smoke'].append(median_smoke)
            results_tap_fast_d['fast_median_sober'].append(median_sober)
            results_tap_fast_d['fast_sd_smoke'].append(sd_smoke)
            results_tap_fast_d['fast_sd_sober'].append(sd_sober)
            results_tap_fast_d['fast_iqr_smoke'].append(iqr_smoke)
            results_tap_fast_d['fast_iqr_sober'].append(iqr_sober)

            tap_Statistic, tap_pvalue, tap_r = self.run_wilcoxon_test(tap_smoke[feature], tap_sober[feature])
            fast_Statistic, fast_pvalue, fast_r = self.run_wilcoxon_test(fast_smoke[feature], fast_sober[feature])
            smoke_Statistic, smoke_pvalue, smoke_Z, smoke_r = self.run_mann_whitney_test(tap_smoke[feature], fast_smoke[feature])
            sober_Statistic, sober_pvalue, sober_Z, sober_r = self.run_mann_whitney_test(tap_sober[feature], fast_sober[feature])

            results_tap_fast['Feature'].append(feature)
            results_tap_fast['tap_Statistic'].append(tap_Statistic)
            results_tap_fast['tap_p-value'].append(tap_pvalue)
            results_tap_fast['tap Effect size'].append(tap_r)
            results_tap_fast['fast_Statistic'].append(fast_Statistic)
            results_tap_fast['fast_p-value'].append(fast_pvalue)
            results_tap_fast['fast Effect size'].append(fast_r)
            results_tap_fast['smoke_Statistic'].append(smoke_Statistic)
            results_tap_fast['smoke_p-value'].append(smoke_pvalue)
            results_tap_fast['smoke_Z'].append(smoke_Z)
            results_tap_fast['smoke Effect size'].append(smoke_r)
            results_tap_fast['sober_Statistic'].append(sober_Statistic)
            results_tap_fast['sober_p-value'].append(sober_pvalue)
            results_tap_fast['sober_Z'].append(sober_Z)
            results_tap_fast['sober Effect size'].append(sober_r)

        reject, pvals_corrected, _, _ = multipletests(results_tap_fast['tap_p-value'], method='holm', alpha=self.alpha)
        results_tap_fast['tap_p-value_corrected'] = pvals_corrected
        results_tap_fast['reject_tap'] = reject

        reject, pvals_corrected, _, _ = multipletests(results_tap_fast['fast_p-value'], method='holm', alpha=self.alpha)
        results_tap_fast['fast_p-value_corrected'] = pvals_corrected
        results_tap_fast['reject_fast'] = reject

        reject, pvals_corrected, _, _ = multipletests(results_tap_fast['smoke_p-value'], method='holm', alpha=self.alpha)
        results_tap_fast['smoke_p-value_corrected'] = pvals_corrected
        results_tap_fast['reject_smoke'] = reject

        reject, pvals_corrected, _, _ = multipletests(results_tap_fast['sober_p-value'], method='holm', alpha=self.alpha)
        results_tap_fast['sober_p-value_corrected'] = pvals_corrected
        results_tap_fast['reject_sober'] = reject

        results_tap_fast = pd.DataFrame(results_tap_fast)
        results_tap_fast.to_excel(os.path.join('Derivatives/all_subs/SemNetAnalysis/group_analyses_results', 'TAP_vs_FAST_group_level.xlsx'), index=False)
        results_tap_fast_d = pd.DataFrame(results_tap_fast_d)
        results_tap_fast_d.to_excel(os.path.join('Derivatives/all_subs/SemNetAnalysis/group_analyses_results', 'TAP_vs_FAST_descriptive_stats.xlsx'), index=False)


    # Function to analyze the data for ses1 and ses2 groups
    def run_ses1_ses2_analysis(self):
        results_ses = {'Feature': [], 'ses1_Statistic': [], 'ses1_p-value': [], 'ses1 Effect size': [], 
                        'ses2_Statistic': [], 'ses2_p-value': [], 'ses2 Effect size': [], 
                        'smoke_Statistic': [], 'smoke_p-value': [], 'smoke_Z': [], 'smoke Effect size': [], 
                        'sober_Statistic': [], 'sober_p-value': [], 'sober_Z': [], 'sober Effect size': []}
        
        results_ses_d = {'Feature': [], 'ses1_mean_smoke':[], 'ses2_mean_sober':[], 'ses2_mean_smoke':[], 'ses1_mean_sober':[],
                         'ses1_median_smoke':[], 'ses2_median_sober':[], 'ses2_median_smoke':[], 'ses1_median_sober':[],
                         'ses1_sd_smoke': [], 'ses2_sd_sober': [], 'ses2_sd_smoke': [], 'ses1_sd_sober': [],
                         'ses1_iqr_smoke': [], 'ses2_iqr_sober': [], 'ses2_iqr_smoke': [], 'ses1_iqr_sober': []}

        ses1_smoke = self.data[(self.data['Smoking'] == 1) & (self.data['Session'] == 1)]
        ses2_sober = self.data[(self.data['Smoking'] == 0) & (self.data['Session'] == 2)]
        ses2_smoke = self.data[(self.data['Smoking'] == 1) & (self.data['Session'] == 2)]
        ses1_sober = self.data[(self.data['Smoking'] == 0) & (self.data['Session'] == 1)]

        for feature in self.features:
            mean_ses1_smoke, median_ses1_smoke, sd_ses1_smoke, iqr_ses1_smoke = self.calculate_descriptive_stats(ses1_smoke[feature])
            mean_ses2_sober, median_ses2_sober, sd_ses2_sober, iqr_ses2_sober = self.calculate_descriptive_stats(ses2_sober[feature])

            results_ses_d['Feature'].append(feature)
            results_ses_d['ses1_mean_smoke'].append(mean_ses1_smoke)
            results_ses_d['ses2_mean_sober'].append(mean_ses2_sober)
            results_ses_d['ses1_median_smoke'].append(median_ses1_smoke)
            results_ses_d['ses2_median_sober'].append(median_ses2_sober)
            results_ses_d['ses1_sd_smoke'].append(sd_ses1_smoke)
            results_ses_d['ses2_sd_sober'].append(sd_ses2_sober)
            results_ses_d['ses1_iqr_smoke'].append(iqr_ses1_smoke)
            results_ses_d['ses2_iqr_sober'].append(iqr_ses2_sober)

            mean_smoke, median_smoke, sd_smoke, iqr_smoke = self.calculate_descriptive_stats(ses2_smoke[feature])
            mean_sober, median_sober, sd_sober, iqr_sober = self.calculate_descriptive_stats(ses1_sober[feature])

            results_ses_d['ses2_mean_smoke'].append(mean_smoke)
            results_ses_d['ses1_mean_sober'].append(mean_sober)
            results_ses_d['ses2_median_smoke'].append(median_smoke)
            results_ses_d['ses1_median_sober'].append(median_sober)
            results_ses_d['ses2_sd_smoke'].append(sd_smoke)
            results_ses_d['ses1_sd_sober'].append(sd_sober)
            results_ses_d['ses2_iqr_smoke'].append(iqr_smoke)
            results_ses_d['ses1_iqr_sober'].append(iqr_sober)

            ses1_Statistic, ses1_pvalue, ses1_r = self.run_wilcoxon_test(ses1_smoke[feature], ses2_sober[feature])
            ses2_Statistic, ses2_pvalue, ses2_r = self.run_wilcoxon_test(ses2_smoke[feature], ses1_sober[feature])
            smoke_Statistic, smoke_pvalue, smoke_Z, smoke_r = self.run_mann_whitney_test(ses1_smoke[feature], ses2_smoke[feature])
            sober_Statistic, sober_pvalue, sober_Z, sober_r = self.run_mann_whitney_test(ses1_sober[feature], ses2_sober[feature])

            results_ses['Feature'].append(feature)
            results_ses['ses1_Statistic'].append(ses1_Statistic)
            results_ses['ses1_p-value'].append(ses1_pvalue)
            results_ses['ses1 Effect size'].append(ses1_r)
            results_ses['ses2_Statistic'].append(ses2_Statistic)
            results_ses['ses2_p-value'].append(ses2_pvalue)
            results_ses['ses2 Effect size'].append(ses2_r)
            results_ses['smoke_Statistic'].append(smoke_Statistic)
            results_ses['smoke_p-value'].append(smoke_pvalue)
            results_ses['smoke_Z'].append(smoke_Z)
            results_ses['smoke Effect size'].append(smoke_r)
            results_ses['sober_Statistic'].append(sober_Statistic)
            results_ses['sober_p-value'].append(sober_pvalue)
            results_ses['sober_Z'].append(sober_Z)
            results_ses['sober Effect size'].append(sober_r)

        reject, pvals_corrected, _, _ = multipletests(results_ses['ses1_p-value'], method='holm', alpha=self.alpha)
        results_ses['ses1_p-value_corrected'] = pvals_corrected
        results_ses['reject_ses1'] = reject

        reject, pvals_corrected, _, _ = multipletests(results_ses['ses2_p-value'], method='holm', alpha=self.alpha)
        results_ses['ses2_p-value_corrected'] = pvals_corrected
        results_ses['reject_ses2'] = reject

        reject, pvals_corrected, _, _ = multipletests(results_ses['smoke_p-value'], method='holm', alpha=self.alpha)
        results_ses['smoke_p-value_corrected'] = pvals_corrected
        results_ses['reject_smoke'] = reject

        reject, pvals_corrected, _, _ = multipletests(results_ses['sober_p-value'], method='holm', alpha=self.alpha)
        results_ses['sober_p-value_corrected'] = pvals_corrected
        results_ses['reject_sober'] = reject

        results_ses_d = pd.DataFrame(results_ses_d)
        results_ses = pd.DataFrame(results_ses)
        results_ses_d.to_excel(os.path.join('Derivatives/all_subs/SemNetAnalysis/group_analyses_results', 'ses1_vs_ses2_descriptive_stats.xlsx'), index=False)
        results_ses.to_excel(os.path.join('Derivatives/all_subs/SemNetAnalysis/group_analyses_results', 'ses1_vs_ses2_group_level.xlsx'), index=False)

    def run_smokeses1_soberses1_analysis(self):
        results_smoke_sober = {'Feature': [], 'Statistic': [], 'p-value': [], 'Effect size': [], 
                               'mean_smoke':[], 'mean_sober':[], 'median_smoke':[], 'median_sober':[],
                                'sd_smoke': [], 'sd_sober': [], 'iqr_smoke': [], 'iqr_sober': []}

        smoke = self.data[(self.data['Smoking'] == 1) & (self.data['Session'] == 2)]
        sober = self.data[(self.data['Smoking'] == 0) & (self.data['Session'] == 2)]

        for feature in self.features:
            statistic, p, z, r = self.run_mann_whitney_test(smoke[feature], sober[feature])
            mean_smoke, median_smoke, sd_smoke, iqr_smoke = self.calculate_descriptive_stats(smoke[feature])
            mean_sober, median_sober, sd_sober, iqr_sober = self.calculate_descriptive_stats(sober[feature])

            results_smoke_sober['Feature'].append(feature)
            results_smoke_sober['Statistic'].append(statistic)
            results_smoke_sober['p-value'].append(p)
            results_smoke_sober['Effect size'].append(r)
            results_smoke_sober['mean_smoke'].append(mean_smoke)
            results_smoke_sober['median_smoke'].append(median_smoke)
            results_smoke_sober['sd_smoke'].append(sd_smoke)
            results_smoke_sober['iqr_smoke'].append(iqr_smoke)
            results_smoke_sober['mean_sober'].append(mean_sober)
            results_smoke_sober['median_sober'].append(median_sober)
            results_smoke_sober['sd_sober'].append(sd_sober)
            results_smoke_sober['iqr_sober'].append(iqr_sober)

        reject, pvals_corrected, _, _ = multipletests(results_smoke_sober['p-value'], method='holm', alpha=self.alpha)
        results_smoke_sober['p-value_corrected'] = pvals_corrected
        results_smoke_sober['reject'] = reject

        results_smoke_sober = pd.DataFrame(results_smoke_sober)
        results_smoke_sober.to_excel(os.path.join('Derivatives/all_subs/SemNetAnalysis/group_analyses_results', 'ses2Smoking_vs_ses2Sober_group_level.xlsx'), index=False)

    

    def run(self):
        # self.run_smoke_sober_analysis()
        # self.run_TAP_FAST_analysis()
        # self.run_ses1_ses2_analysis()
        self.run_smokeses1_soberses1_analysis()

if __name__ == "__main__":
    features = ['Nodes', 'Edges', 'degree_centrality', 'betweenness_centrality', 
                'clustering_coefficient', 'average_path_length',
                'modularity', 'communities', 'word_count']

    # Remove outliers
    analysis_df = pd.read_excel('Derivatives/all_subs/SemNetAnalysis/network_stats.xlsx')
    # Run statistical analysis on the merged data
    analysis = Analysis(analysis_df, features)
    analysis.run() 