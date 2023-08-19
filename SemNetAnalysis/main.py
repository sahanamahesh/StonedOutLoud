"""
File: main.py
Author: Sahana Mahesh
Co-Author: Jen Burrell
Email: sahana.sm61@gmail.com
Date Created: July 20th, 2023
Last Modified: August 2nd, 2023
Description: This script is part of a COGS 402 project (UBC 2023S tearm 1 & 2) conducted at the Cognitive Neuroscience of Thought Laboratory 
             under the supervision of Dr. Kalina Christoff (PI) and Jen Burrell. This script is used to process, analyze, and generate   
             semantic network statistics of participant data. 
"""
from preprocess import Preprocessor
from cooccurrence_network import CooccurrenceNetworkBuilder
from network_stats import NetworkStats
from analysis import Analysis
import os
import networkx as nx
import glob
import pandas as pd
import openpyxl
import msoffcrypto
import io

def main(WDIR, sub_list):

    stats_data_frame = pd.DataFrame()

    # Retrieve a list of paths to all TAP transcriptions
    file_list = sorted(glob.glob(WDIR+"/*/*/Transcriptions/*-TAP.txt"))

    # dictionary to keep track of the sessions for each participant
    participant_sessions = {}

    passwd = 'SOL_2023'

    decrypted_workbook = io.BytesIO()
    with open('Data/Participant_Log.xlsx', 'rb') as file:
        office_file = msoffcrypto.OfficeFile(file)
        office_file.load_key(password=passwd)
        office_file.decrypt(decrypted_workbook)

    data_log = pd.read_excel(decrypted_workbook)

    df1 = data_log.copy()

    #stats_data_frame['Subject'] = (stats_data_frame['Subject'].astype(str)).str.upper()
    df1['Subject'] = (df1['Subject'].astype(str)).str.upper()

    df1['Smoking'] = df1['Smoke Session'].apply(lambda x: 1 if x == 1 else 0)
    df1['TAP'] = df1['Session start task'].apply(lambda x: 1 if x == 'TAP' else 0)

    #stats_data_frame = pd.merge(stats_data_frame, df1[['Subject', 'Smoking', 'TAP']], on=['Subject'], how='left')
    
    # Iterate over each file in file_list
    for file in file_list:

        # Skip files that contain '.log' in filename
        if '.log' in file:
            continue

        participant_id = file.split('/')[-4]
        session = file.split('/')[-3]

        if participant_id not in participant_sessions:
            participant_sessions[participant_id] = [session]
        else:
            participant_sessions[participant_id].append(session)


    for participant_id, sessions in participant_sessions.items():

        if 'ses-1' not in sessions or 'ses-2' not in sessions:
            continue

        participant_files = [file for file in file_list if participant_id in file]

        for file in participant_files:

            session = file.split('/')[-3]
            sub_num = participant_id[4:]
            ses_num = session[-1]

            # Check if it's a smoking session for the subject
            subject_smoking_values = df1.loc[df1['Subject'] == sub_num, 'Smoking']
            smoking_value = subject_smoking_values.values[0] if not subject_smoking_values.empty else 0
            #smoking = df1['Smoking'].iloc[0]

            # Seperate subject filename from path and store it in sub
            sub = file.rsplit('/',1)[1].split('_',1)[0]

            # Skip subjects (sub) that are in sub_list
            if sub.split('-')[-1] in sub_list:
                continue

            # file_out holds the name for the cleaned text file    
            file_out = file.split('.')[0] + "_clean.txt" 

            # output_dir holds the output directory name
            output_dir = file.split('Transcriptions')[0] + "Networks" 

            # Clean text
            preprocessor = Preprocessor(file)
            word_count = preprocessor.get_word_count()
            tokens = preprocessor.preprocess()

            # Add cleaned text to file_out
            preprocessor.write_to_file(file_out)  

            # Estimate network
            network_builder = CooccurrenceNetworkBuilder(tokens)
            network = network_builder.build_cooccurrence_network()

            # Group level network statistics
            network_stats = NetworkStats(network)
            stats = network_stats.calculate_stats()
            
            stats['word_count'] = word_count
            stats['Subject'] = sub_num 
            stats['Session'] = ses_num 
            
        
            if not os.path.exists('Derivatives/all_subs/SemNetAnalysis/sober_networks'):
                os.makedirs('Derivatives/all_subs/SemNetAnalysis/sober_networks')

            if not os.path.exists('Derivatives/all_subs/SemNetAnalysis/cannabis_networks'):
                os.makedirs('Derivatives/all_subs/SemNetAnalysis/cannabis_networks')

            if not os.path.exists('Derivatives/all_subs/SemNetAnalysis/group_analyses_results'):
                os.makedirs('Derivatives/all_subs/SemNetAnalysis/group_analyses_results')
            
            # Subject level network statistics
            # individual_stats = network_stats.calculate_individual_stats()
            # individual_stats_df = pd.DataFrame(individual_stats)
            # individual_stats_df['Subject'] = sub_num
            # individual_stats_df['Session'] = ses_num

            #
            # Add the 'Smoking' information to the filename
            if smoking_value == 1 and ses_num == '1':
                #individual_stats_filename = f'sub-{sub_num}_ses-{ses_num}_smoking.xlsx'
                individual_network_filename = f'sub-{sub_num}_ses-{ses_num}_smoking.gexf'
                stats['Smoking'] = 1
                #individual_stats_df.to_excel(os.path.join('stats/subject_stats/cannabis', individual_stats_filename), index = True)
                nx.write_gexf(network, os.path.join('Derivatives/all_subs/SemNetAnalysis/cannabis_networks', individual_network_filename))
            elif smoking_value == 0 and ses_num == '2':
                #individual_stats_filename = f'sub-{sub_num}_ses-{ses_num}_smoking.xlsx'
                individual_network_filename = f'sub-{sub_num}_ses-{ses_num}_smoking.gexf' #f'{sub_num}_smoking-ses{ses_num}.xlsx'
                stats['Smoking'] = 1
                #individual_stats_df.to_excel(os.path.join('stats/subject_stats/cannabis', individual_stats_filename), index = True)
                nx.write_gexf(network, os.path.join('Derivatives/all_subs/SemNetAnalysis/cannabis_networks', individual_network_filename))
            else:
                individual_network_filename = f'sub-{sub_num}_ses-{ses_num}_smoking.gexf'
                #individual_stats_filename = f'sub-{sub_num}_ses-{ses_num}_sober.xlsx'
                stats['Smoking'] = 0
                #individual_stats_df.to_excel(os.path.join('stats/subject_stats/sober', individual_stats_filename), index = True)
                nx.write_gexf(network, os.path.join('Derivatives/all_subs/SemNetAnalysis/sober_networks', individual_network_filename))

            stats_data_frame = pd.concat([stats_data_frame, stats], ignore_index = True)
            stats_data_frame.to_excel(os.path.join('Derivatives/all_subs/SemNetAnalysis', 'network_stats.xlsx'), index = False)

            # individual_stats_df.to_excel(os.path.join('stats/subject_stats', individual_stats_filename), index = True)

            # Create a new folder to store the estimated network if it doesnt exist
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            # Save the network as a .gexf file in the Network folder
            nx.write_gexf(network, os.path.join(output_dir, file.replace('.txt', '-Network.gexf').rsplit('/', 1)[1]))


    stats_data_frame['Subject'] = (stats_data_frame['Subject'].astype(str)).str.upper()
    df1['Subject'] = (df1['Subject'].astype(str)).str.upper() 

    # df1['Smoking'] = df1['Smoke Session'].apply(lambda x: 1 if x == 1 else 2)
    df1['TAP'] = df1['Session start task'].apply(lambda x: 1 if x == 'TAP' else 0)

    stats_data_frame = pd.merge(stats_data_frame, df1[['Subject', 'TAP']], 
                            on=['Subject'], how='left')
    
    stats_data_frame.to_excel(os.path.join('Derivatives/all_subs/SemNetAnalysis', 'network_stats.xlsx'), index = False)

    features = ['Nodes', 'Edges', 'degree_centrality', 'betweenness_centrality', 
                'clustering_coefficient', 'average_path_length',
                'modularity', 'communities', 'word_count']

    # Remove outliers
    analysis_df = pd.read_excel('Derivatives/all_subs/SemNetAnalysis/network_stats.xlsx')
    # Run statistical analysis on the merged data
    analysis = Analysis(analysis_df, features)
    analysis.run() # This saves the analysis results to new Excel files

if __name__ == "__main__":
    WDIR = '/Users/sahanamahesh/Library/CloudStorage/OneDrive-SharedLibraries-UBC/Burrell, Jen - StonedOutLoud/Derivatives'
    sub_list = ['P001', 'P003', 'P007', 'P004', '007', '011', '017', '019', '028', '029', '045']
    main(WDIR, sub_list)