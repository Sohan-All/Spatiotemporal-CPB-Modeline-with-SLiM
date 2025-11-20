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

import warnings

#WARNING: don't run this file in VSCode. Run it in the terminal instead.

def main():
    #hehehe
    warnings.filterwarnings("ignore")
    
    #Ask for input on number of clusters
    num_clusters = int(input("Enter the number of clusters (default 99): ").strip() or 99)
    
    #Query for a migration rates modifier
    migration_rates_modifier = int(input("Enter the migration rates modifier (default 1000): ").strip() or 1000)
    
    #Query for mutation rate
    mutation_rate = float(input("Enter the mutation rate (default 1e-7): ").strip() or 1e-7)

    #Query for recombination rate
    recombination_rate = float(input("Enter the recombination rate (default 1e-8): ").strip() or 1e-8)
    
    #Query for population modifier
    population_modifier = float(input("Enter the total population size (default 10000): ").strip() or 10000)
    
    #Start by reading the data from final_data_for_modeling.csv
    print("Setting up data for simulations...")
    field_data = CollectData.read_csv('../data/final_data_for_modeling.csv')
    
    #Cluster the coordinates using KMeans
    GenerateClusterData.cluster_coordinates(field_data, n_clusters=num_clusters, iters=2000, random_state=random.randint(0, 1000))
    
    #Put the data for clusters into a list of Cluster objects
    clusters = GenerateClusterData.populate_cluster_objects(field_data, estimate_data=True)
    
    #Assign genomes to clusters based on the specifier matrix for a certain year
    GenerateClusterData.assign_genomes_to_clusters(clusters)
    
    #Generate a distance matrix for the clusters
    distances = GenerateClusterData.create_cluster_distance_matrix(clusters, output_path='../data/cluster_distances.csv')   
     
    #Save the cluster data to a CSV file
    GenerateClusterData.cluster_data_to_csv(clusters, output_path='../data/cluster_data.csv')
        
    #Generate migration rates based on the cluster distance matrix
    GenerateSimulationParams.determine_migration_rates(distances, modifier=migration_rates_modifier, output_path='../data/migration_rates.csv')
    
    #Run the SLiM simulation to create the tree sequence
    print("Running SLiM simulation...")
    subprocess.run(['slim', '-l', '0', '-d', f'POPMULT={population_modifier}', '../SLiM Code/CPBSampleSim.slim'])
    
    #Does recapitation and mutation addition, then gets diversity and divergence statistics
    print("Recapitating tree sequence...")
    AnalyzeTreeSeq.analyze_tree_sequence(mutation_rate=mutation_rate, recombination_rate=recombination_rate)
    print("Successfully generated diversity and divergence statistics from tree sequence.")
    
    
    

    
main()
