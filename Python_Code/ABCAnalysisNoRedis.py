import csv
import numpy as np

from pathlib import Path
from scipy import stats

import Main

# Define prior distributions using scipy.stats
prior_distributions = {
    "m": stats.lognorm(s=1.5, scale=np.exp(np.log(0.0001))),
    "pop": stats.expon(loc=2000, scale=50000),
    "numClusters": stats.randint(1, 4)  # randint(1, 4) gives 1, 2, or 3
}



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
    pop = int(np.floor(parameter["pop"]))
    numClusters = parameter["numClusters"] * 33  #scale to 33, 66, or 99
    
    #Run the model TODO:change silent to true for actual runs
    Main.main(num_clusters=numClusters, migration_rates_modifier=m, population_modifier=pop, silent=False)
    
    
    #Read in the output data
    outDict = {}
    
    for year in ["2015", "2019", "2023"]:
        with open(Path(f"../data/Output_data/diversities_{year}.csv"), mode='r', newline='', encoding='utf-8') as csvfile:
            div2015 = csv.DictReader(csvfile)
            diversities_list = []
            for row in div2015:
                value = next(iter(row.values()))
                if value is not None and value.strip() != "":
                    diversities_list.append(float(value.strip()))
            diversities = np.array(diversities_list)
            outDict[f"{year}_diversity"] = diversities
        
        with open(Path(f"../data/Output_data/divergences_{year}.csv"), mode='r', newline='', encoding='utf-8') as csvfile:
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
            outDict[f"{year}_divergence"] = divergences
        

    return outDict



def calculate_losses(x, x0):
    '''
    Calculate loss metrics comparing observed and simulated data.
    x is observed, x0 is simulated.
    
    Returns a dictionary containing:
    - diversity_loss_{year} for each year
    - divergence_loss_{year} for each year
    - total_loss
    '''
    losses = {}
    total_loss = 0
    
    for year in ["2015", "2019", "2023"]:
        diversity_key = f"{year}_diversity"
        divergence_key = f"{year}_divergence"
        
        # Pi (diversity) distance
        diversity_loss = 0
        for i in range(len(x[diversity_key])):
            pi_distance = abs(x[diversity_key][i] - x0[diversity_key][i])
            pi_distance *= len(x[diversity_key])  # weight pi distance equally to fst distance
            diversity_loss += pi_distance
        losses[f"diversity_loss_{year}"] = diversity_loss
        total_loss += diversity_loss
        
        # Fst (divergence) distance
        divergence_loss = 0
        for i in range(len(x[divergence_key])):
            for k in range(len(x[divergence_key])):
                if i != k:
                    fst_distance = abs(x[divergence_key][i][k] - x0[divergence_key][i][k])
                    divergence_loss += fst_distance
        losses[f"divergence_loss_{year}"] = divergence_loss
        total_loss += divergence_loss
    
    losses["total_loss"] = total_loss
    return losses

def getObservedData():
    outDict = {}
    
    for year in ["2015", "2019", "2023"]:
        path = f"../data/empiricalStats/averaged_pi_{year}.csv"
        with open(Path(path), mode='r', newline='', encoding='utf-8') as csvfile:
            div2015 = csv.DictReader(csvfile)
            diversities_list = []
            for row in div2015:
                value = next(iter(row.values()))
                if value is not None and value.strip() != "":
                    diversities_list.append(float(value.strip()))
            diversities = np.array(diversities_list)
            outDict[f"{year}_diversity"] = diversities
        
        path = f"../data/empiricalStats/averaged_dxy_{year}.csv"
        with open(Path(path), mode='r', newline='', encoding='utf-8') as csvfile:
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
            outDict[f"{year}_divergence"] = divergences
            
    return outDict
    

def sample_prior():
    '''
    Sample parameters from the prior distributions.
    '''
    return {
        "m": prior_distributions["m"].rvs(),
        "pop": prior_distributions["pop"].rvs(),
        "numClusters": prior_distributions["numClusters"].rvs()
    }


def run_abc_simulation(num_iterations, output_csv="../out/abc_results.csv"):
    '''
    Run ABC simulations by repeatedly sampling from the prior and computing losses.
    Results are appended to a CSV file.
    
    :param num_iterations: Number of iterations to run
    :param output_csv: Path to the output CSV file
    '''

    observed_data = getObservedData()
    
    # Determine if we need to write the header
    csv_exists = Path(output_csv).exists()
    
    # Define CSV columns
    fieldnames = ["iteration", "m", "pop", "numClusters",
                  "diversity_loss_2015", "diversity_loss_2019", "diversity_loss_2023",
                  "divergence_loss_2015", "divergence_loss_2019", "divergence_loss_2023",
                  "total_loss"]
    
    with open(output_csv, mode='a', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        # Write header only if file is new
        if not csv_exists:
            writer.writeheader()
        
        for iteration in range(num_iterations):
            parameters = sample_prior()
            print(
                f"Running iteration {iteration + 1}/{num_iterations} "
                f"with m={parameters['m']:.6g}, pop={int(np.floor(parameters['pop']))}, "
                f"numClusters={parameters['numClusters'] * 33}..."
            )
            
            # Sample from prior
            parameters = sample_prior()
            
            try:
                # Run the model
                simulated_data = model(parameters)
                
                # Calculate losses
                losses = calculate_losses(observed_data, simulated_data)
                
                # Prepare row for CSV
                row = {
                    "iteration": iteration,
                    "m": parameters["m"],
                    "pop": parameters["pop"],
                    "numClusters": parameters["numClusters"],
                    "diversity_loss_2015": losses["diversity_loss_2015"],
                    "diversity_loss_2019": losses["diversity_loss_2019"],
                    "diversity_loss_2023": losses["diversity_loss_2023"],
                    "divergence_loss_2015": losses["divergence_loss_2015"],
                    "divergence_loss_2019": losses["divergence_loss_2019"],
                    "divergence_loss_2023": losses["divergence_loss_2023"],
                    "total_loss": losses["total_loss"]
                }
                
                # Append to CSV
                writer.writerow(row)
                csvfile.flush()  # Ensure data is written immediately
                
                print(f"  Total loss: {losses['total_loss']:.6f}")
                
            except Exception as e:
                print(f"  Error in iteration {iteration}: {e}")
                continue
    
    print(f"Simulation {iteration + 1} complete. Results saved to {output_csv}")


if __name__ == "__main__":
    # Run 100 iterations of ABC sampling
    # Modify num_iterations as desired
    run_abc_simulation(num_iterations=1, output_csv="../out/abc_results.csv")