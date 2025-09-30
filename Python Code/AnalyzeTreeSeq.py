import tskit
import pandas as pd
import math


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
#ts = ts.simplify(samples=genome_indicies)



# Calculate diversity for each population in genome_indicies


diversities = []
for idx in genome_indicies:
    # Get the sample node for this population
    pop_samples = ts.samples(population=idx)
    #calculate diversity and append to list
    pi = ts.diversity([pop_samples])
    diversities.append(pi)
    
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
            print(p1)
            print(p2)
            print(ts.num_sites)
            fst = ts.Fst([p1, p2])
            divergences[i].append(fst)

for divergence in divergences:
    print(divergence)