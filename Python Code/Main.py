import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt

import DataWrappers
import CollectData
import GenerateSimulationParams
import GenerateClusterData
import AnalyzeTreeSeq

import subprocess

#WARNING: don't run this file in VSCode. Run it in the terminal instead.

def main():
    #Start by reading the data from final_data_for_modeling.csv
    field_data = CollectData.read_csv('../data/final_data_for_modeling.csv')
    
    #Ask for input on number of clusters
    num_clusters = int(input("Enter the number of clusters (default 99): ").strip() or 99)
    
    #Cluster the coordinates using KMeans
    GenerateClusterData.cluster_coordinates(field_data, n_clusters=num_clusters, iters=2000, random_state=random.randint(0, 1000))
    
    #Put the data for clusters into a list of Cluster objects
    clusters = GenerateClusterData.populate_cluster_objects(field_data, estimate_data=True)
    
    #Assign genomes to clusters based on the specifier matrix for a certain year
    GenerateClusterData.assign_genomes_to_clusters(clusters, specifier_matrix="../data/Genetic_Data/specifier_matrix_2023.csv")
    
    #TODO: fix logic to do with choosing the genomes for clusters and make sure that all years are done at the same time.

    #Generate a distance matrix for the clusters
    distances = GenerateClusterData.create_cluster_distance_matrix(clusters, output_path='../data/cluster_distances.csv')   
     
    #Save the cluster data to a CSV file
    GenerateClusterData.cluster_data_to_csv(clusters, output_path='../data/cluster_data.csv')
    
    #Query for a migration rates modifier
    migration_rates_modifier = int(input("Enter a migration rates modifier (default 10000): ").strip() or 10000)
    
    #Generate migration rates based on the cluster distance matrix
    migration_rates = GenerateSimulationParams.determine_migration_rates(distances, modifier=migration_rates_modifier, output_path='../data/migration_rates.csv')
    
    #Run the SLiM simulation to create the tree sequence
    subprocess.run(['slim', '../SLiM Code/CPBSampleSim.slim'])
    
    #Does recapitation and mutation addition, then gets diversity and divergence statistics
    AnalyzeTreeSeq.analyze_tree_sequence()
    
    
    

    
main()
