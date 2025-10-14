import tskit
import msprime, pyslim
import pandas as pd
import math
import numpy as np
import csv

def analyze_tree_sequence():
    # Load the cluster_data CSV file
    cluster_data = pd.read_csv("../data/cluster_data.csv")

    assignments = cluster_data['Genome Assignment']

    total_assignments = int(max(assignments.dropna())) 

    genome_indicies = [-1] * (total_assignments+1)

    for i in range(len(assignments)):
        if math.isnan(assignments[i]) == False:
            index = int(assignments[i])
            if genome_indicies[index] == -1:
                genome_indicies[index] = i

    # Load the tree sequence file
    ts = tskit.load("../out/simTreeSeq.trees")


    diversities = []
    for idx in genome_indicies:
        # Get the sample node for this population
        pop_samples = ts.samples(population=idx)
        #calculate diversity and append to list
        pi = ts.diversity([pop_samples])
        diversities.append(pi[0])
    diversities = np.array(diversities)
        
    print("Diversity (pi) for each genome assignment before recapitation and mutation addition:")
    print(diversities)


    ts = pyslim.recapitate(ts, recombination_rate=1e-8, ancestral_Ne=6700)
    #ts = ts.simplify(samples=genome_indicies)
    next_id = pyslim.next_slim_mutation_id(ts)
    ts = msprime.sim_mutations(
            ts,
            rate=1e-8,
            model=msprime.SLiMMutationModel(type=0, next_id=next_id),
            keep=True,
    )

    print(f"The tree sequence now has {ts.num_mutations} mutations,\n"
        f"and mean pairwise nucleotide diversity is {ts.diversity():0.3e}.")


    # Calculate diversity for each population in genome_indicies


    diversities = []
    for idx in genome_indicies:
        # Get the sample node for this population
        pop_samples = ts.samples(population=idx)
        #calculate diversity and append to list
        pi = ts.diversity([pop_samples])
        diversities.append(pi[0])
    diversities = np.array(diversities)
        
    print("Diversity (pi) for each genome assignment:")
    print(diversities)
        
    divergences = []

    for i in range(len(genome_indicies)):
        divergences.append([])
        for j in range(len(genome_indicies)):
            if i == j:
                divergences[i].append(0)
            #elif i > j:
            #    divergences[i].append(divergences[j][i])
            else:
                p1 = ts.samples(population=genome_indicies[i])
                p2 = ts.samples(population=genome_indicies[j])
                divergence = ts.divergence([p1, p2])
                divergences[i].append(divergence.item())

    # for divergence in divergences:
    #     print(divergence)
    
    
    #write the data
    with open("../data/diversitiescsv", "w", newline="") as f:
        writer = csv.writer(f)
        for item in diversities:
            writer.writerow([item])
        
    with open("../data/divergences.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(divergences)
        
    
        
#analyze_tree_sequence()