import pandas as pd
from collections import Counter
import csv

CUTOFF = 10000

def create_specifier_matrix(year, genetic_data="../data/Genetic_Data/CPB_genetic_metadata.csv", output_path="../data/Genetic_Data/specifier_matrix.csv"):
    """
    Create a specifier matrix for the given year from genetic data. 
    This matrix will have rows with population names, latitude, longitude, and indices of genetic data.
    
    Parameters:
    year (int): The year for which to create the specifier matrix.
    genetic_data (str): Path to the genetic data CSV file.
    output_path (str): Path to save the specifier matrix CSV file.
    
    Returns:
    None
    """
    # Read the genetic data
    df = pd.read_csv(genetic_data)
    
    matrix = []
    currPop = None
    
    for index, row in df.iterrows():
        # You can access each row's data using 'row'
        if row['Year'] != year:
            continue
        # print(row)
        # print("----------------------------")
        if currPop is None or currPop != row['Population']:
            matrix.append([row['Population'], row['Latitude'], row['Longitude']])
            currPop = row['Population']
        else:
            matrix[-1].append(index)
        
        #print(matrix)
        # # Convert matrix to DataFrame and save to CSV
        specifier_df = pd.DataFrame(matrix)
        specifier_df.to_csv(output_path, index=False, header=False)

    

def generate_Ref_FASTA_from_genolike(cutoff, filePath="../../../Data for modeling/cpbWGS_genolike_chr9.cpbWGS_genolike_chr9.beagle.gz.phased", output_path_fasta="../data/Genetic_Data/refFull.fasta", output_path_ids="../data/Genetic_Data/FastaIDs.csv"):
    """
    Generate a reference FASTA file from genetic data.
    This function reads the genetic data and creates a FASTA file.
    
    Parameters:
    cutoff (int): The cutoff value for filtering genetic data. If the percentage of similar base pairs is above this cutoff,
            it is not included in the FASTA file.
    filePath (str): Path to the genetic data file. This file is not in the github repository and should be provided separately.
    
    Returns:
    None
    """
    with open(filePath, 'r') as f:
        #lines = f.readlines()
        lines = [next(f) for _ in range(CUTOFF)] #TODO: remove this line and uncomment the previous line to read the full file
        
    #print("read done")

    matrix = [line.strip().split() for line in lines if line.strip()]
    
    # Remove the first row and first column from the matrix which just stores headers
    matrix = [row[1:] for row in matrix[1:]]
    fastaStr = ""
    
    # Extract the first column into a separate list
    ids = [row[0] for row in matrix]
    # Remove the first column from each row in the matrix
    matrix = [row[1:] for row in matrix]
    
    usedIDs = []
    
    #print("started creating fasta")
    
    # Add the most common base from each row to the FASTA string
    line_length = 80
    for i in range(len(matrix)):
        most_common, count = Counter(matrix[i]).most_common(1)[0]
        if count / len(matrix[i]) <= cutoff:
            fastaStr += most_common
            usedIDs.append(ids[i])
            if len(fastaStr) % line_length == 0:
                fastaStr += "\n"
    
    #print("started writing fasta")
    
    # Convert numeric bases to nucleotides
    base_map = {'0': 'A', '1': 'C', '2': 'G', '3': 'T'}
    fastaStr = ''.join([base_map.get(base, base) for base in fastaStr])
    
    with open(output_path_fasta, 'w') as fasta_file:
        fasta_file.write(">refFull\n")
        fasta_file.write(fastaStr + "\n")
    
    # Write usedIDs to a CSV file with one column
    with open(output_path_ids, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID'])
        for uid in usedIDs:
            writer.writerow([uid])
        
    
#Create a method here that generates mutation lists for each individual in the relevant year
#The method should take in a year, use the refFull.fasta file, and the specifier matrix to generate multiple lists of mutations
#for each relevant individual that was sequenced in that year.

def generateMutationLists(year, specifier_matrix="../data/Genetic_Data/specifier_matrix_2023.csv", 
                          genetic_data="../../../Data for modeling/cpbWGS_genolike_chr9.cpbWGS_genolike_chr9.beagle.gz.phased",
                          ref_fasta="../data/Genetic_Data/refFull.fasta", 
                          fasta_ids="../data/Genetic_Data/FastaIDs.csv",
                          output_path="../data/Genetic_Data/mutation_lists.csv"):
    '''
    Generate mutation lists for each group of individuals in the same location in the specified year.
    This function reads the reference FASTA file, the specifier matrix, and the genetic data file,
    and creates a CSV file with mutation lists for each location.
    
    Parameters:
    year (int): The year for which to generate mutation lists.
    specifier_matrix (str): Path to the specifier matrix CSV file.
    genetic_data (str): Path to the genetic data file.
    ref_fasta (str): Path to the reference FASTA file.
    fasta_ids (str): Path to the FastaIDs CSV file.
    output_path (str): Path to save the mutation lists CSV file.
    Returns:
    None
    '''
    
    # Read the reference FASTA file and extract the base pair string
    with open(ref_fasta, 'r') as f:
        lines = f.readlines()
        # Remove the first line (comment) and join the rest into a single string without newlines
        ref_sequence = ''.join(line.strip() for line in lines[1:])
        
    # Read the FastaIDs file and put each ID into a list
    fasta_id_list = []
    with open(fasta_ids, 'r') as f:
        next(f)  # Skip header
        for line in f:
            fasta_id_list.append(line.strip())
            
    # Read the specifier matrix and extract genome indexes for the relevant year
    specifier_df = pd.read_csv(specifier_matrix, header=None)
    # The genome indexes start from the 4th column (index 3)
    genome_indexes_matrix = specifier_df.iloc[:, 3:].values.tolist()
    
    # Read the genetic data file as a matrix
    with open(genetic_data, 'r') as f:
        #genetic_lines = [line.strip().split() for line in f if line.strip()]
        genetic_lines = []
        with open(genetic_data, 'r') as f:
            for i, line in enumerate(f):
                if i >= CUTOFF:
                    break
                if line.strip():
                    genetic_lines.append(line.strip().split())
        

    mutation_lists = [[] for _ in genome_indexes_matrix]

    #print(genome_indexes_matrix)

    # Iterate over each line of the genetic data. 
    # This part condenses all genetic data from one spatiotemporal point into a single list of mutations 
    # based on the most common base at each location. It compares that most common base to the reference sequence
    # and if they differ, it adds the mutation to the corresponding list.
    counter = 0
    for line in genetic_lines[1:]:
        #print(counter)
        if line[1] != fasta_id_list[counter]:
            continue
        for i in range(len(genome_indexes_matrix)):
            bases = []
            for index in range(len(genome_indexes_matrix[i])):
                # Check for NaN and break if found
                if pd.isna(genome_indexes_matrix[i][index]):
                    break
                bases.append(line[int(genome_indexes_matrix[i][index] * 2 + 1)])
                bases.append(line[int(genome_indexes_matrix[i][index] * 2 + 2)])
            most_common_base = Counter(bases).most_common(1)[0][0]
            if int_to_base(int(most_common_base)) != ref_sequence[counter]:
                #print(bases, most_common_base, ref_sequence[counter])
                mutation_lists[i].append((counter, int_to_base(int(most_common_base))))
        counter += 1
        if counter >= len(ref_sequence):
            break
        
        
        
    for i, mutation_list in enumerate(mutation_lists):
        print(f"Mutation list {i+1}:")
        print(mutation_list)
        print('-' * 40)
    
    # Output each element of mutation_lists to a CSV such that each list within the 2D list takes up two columns [mutation idx, base]
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        header = []
        for i in range(len(mutation_lists)):
            header.extend([f"Index_{i+1}", f"Base_{i+1}"])
        writer.writerow(header)
        # Find the max length among all mutation lists
        max_len = max(len(lst) for lst in mutation_lists)
        # Write rows
        for row_idx in range(max_len):
            row = []
            for lst in mutation_lists:
                if row_idx < len(lst):
                    row.extend([lst[row_idx][0], lst[row_idx][1]])
                else:
                    row.extend(['', ''])
            writer.writerow(row)
    
    
def int_to_base(num):
    """
    Convert an integer to a nucleotide base. 0=A, 1=C, 2=G, 3=T.
    
    Parameters:
    num (int): The integer to convert (0, 1, 2, or 3).
    
    Returns:
    str: The corresponding nucleotide base ('A', 'C', 'G', or 'T').
    """
    return {0: 'A', 1: 'C', 2: 'G', 3: 'T'}.get(num, '')

    
    
#create_specifier_matrix(2023, genetic_data="../data/Genetic_Data/CPB_genetic_metadata.csv", output_path="../data/Genetic_Data/specifier_matrix_2023.csv")
generate_Ref_FASTA_from_genolike(0.95)
generateMutationLists(2023)