import csv
import os
import tempfile

import matplotlib.pyplot as plt
import numpy as np

import pyabc

from pyabc.transition import MultivariateNormalTransition

import Main

pyabc.settings.set_figure_params('pyabc')  # for beautified plots

prior = pyabc.Distribution(
    m=pyabc.RV("lognorm", np.log(0.0001), 1.5),
    pop=pyabc.RV("expon", loc=2000, scale=50000),
    numClusters=pyabc.RV("randint", 1, 3)
)




def model(parameter):
    '''
    The model function that runs the SLiM simulation with the given parameters: 
    1. migration rate (m)
    2. population size (pop)
    3. number of clusters (numClusters)   
    
    :param parameter: This is a dictionary containing the parameters for the simulation.
    '''
    
    #Get the parameters
    m = parameter["m"]
    pop = parameter["pop"]
    numClusters = parameter["numClusters"] * 33  #scale to 33, 66, or 99
    
    #Run the model
    Main.main(num_clusters=numClusters, migration_rates_modifier=m, population_modifier=pop, silent=True)
    
    
    #Read in the output data
    outDict = {"2015": {}, "2019": {}, "2023": {}}
    
    for year in ["2015", "2019", "2023"]:
        with open(f"..\\data\\Output_data\\diversities_{year}.csv", mode='r', newline='', encoding='utf-8') as csvfile:
            div2015 = csv.DictReader(csvfile)
            diversities_list = []
            for row in div2015:
                value = next(iter(row.values()))
                if value is not None and value.strip() != "":
                    diversities_list.append(float(value.strip()))
            diversities = np.array(diversities_list)
            outDict[year]["diversities"] = diversities
        
        with open(f"..\\data\\Output_data\\divergences_{year}.csv", mode='r', newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            matrix = []
            for row in reader:
                if not row:
                    continue
                row_vals = []
                for val in row:
                    v = val.strip()
                    if v == "":
                        row_vals.append(np.nan)
                    else:
                        row_vals.append(float(v))
            matrix.append(row_vals)
            divergences = np.array(matrix, dtype=float)
            outDict[year]["divergences"] = divergences
        

    return outDict



def distance(x, x0):
    '''
    Distance function for comparing simulated and observed data.
    x is observed, x0 is simulated.
    Pi is nucleotide diversity given in the form of a list of values for each population.
    Fst is genetic differentiation given in the form of a matrix of values between populations.
    
    x should be a dictionary with keys "2015", "2019", "2023". Each of these keys has values: "diversities" and "divergences"
    which are the pi and fst values in list and matrix form, respectively.
    
    Each of these values is an array of three of these values, one for each year (2015, 2019, 2023).
    
    Returns a single distance value.
    '''
    total_distance = 0
    
    for year in ["2015", "2019", "2023"]:    
        # Pi distance
        for i in range(len(x[year]["diversities"])):
            pi_distance = abs(x[year]["diversities"][i] - x0[year]["diversities"][i])
            pi_distance *= len(x[year]["diversities"]) #weight pi distance equally to fst distance
            total_distance += pi_distance
        # Fst distance
        for i in range(len(x[year]["divergences"])):
            for k in range(len(x[year]["divergences"])):
                if i != k:
                    fst_distance = abs(x[year]["divergences"][i][k] - x0[year]["divergences"][i][k])
                    total_distance += fst_distance
                    
    return total_distance



abc = pyabc.ABCSMC(model, prior, distance, population_size=500, transitions=MultivariateNormalTransition())

db_path = os.path.join(tempfile.gettempdir(), "test.db")
observation = {"singlesChosen": 7}
abc.new("sqlite:///" + db_path, observation)

history = abc.run(minimum_epsilon=2, max_nr_populations=5)
df, weights = history.get_distribution()
df.to_csv("posterior_samples.csv", index=False)

# df, w = history.get_distribution()
# posterior_mean = (df * w[:, None]).sum()
# print("Posterior mean:", posterior_mean)