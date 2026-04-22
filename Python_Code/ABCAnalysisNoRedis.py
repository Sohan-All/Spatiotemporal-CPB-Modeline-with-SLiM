import csv
import numpy as np

from pathlib import Path
from scipy import stats

import Main

# Default parameter values for log-normal distributions
DEFAULT_MUTATION_RATE = 2.1e-9
DEFAULT_RECOMBINATION_RATE = 2.75e-6

# Define prior distributions using scipy.stats
prior_distributions = {
    "m": stats.lognorm(s=1.5, scale=np.exp(np.log(0.0001))),
    "pop": stats.expon(loc=2000, scale=50000),
    "numClusters": stats.randint(1, 4),  # randint(1, 4) gives 1, 2, or 3
    "mutation_rate": stats.lognorm(s=0.5, scale=DEFAULT_MUTATION_RATE),  # log-normal around default
    "recombination_rate": stats.lognorm(s=0.5, scale=DEFAULT_RECOMBINATION_RATE)  # log-normal around default
}



def model(parameter):
    '''
    The model function that runs the SLiM simulation with the given parameters: 
    1. migration rate (m)
    2. population size (pop)
    3. number of clusters (numClusters)
    4. mutation rate (mutation_rate)
    5. recombination rate (recombination_rate)
    
    :param parameter: This is a dictionary containing the parameters for the simulation.
    '''
    
    #Get the parameters
    m = parameter.get("m", prior_distributions["m"].rvs())
    pop = int(np.floor(parameter.get("pop", prior_distributions["pop"].rvs())))
    numClusters = parameter.get("numClusters", prior_distributions["numClusters"].rvs()) * 33  #scale to 33, 66, or 99
    mutation_rate = parameter.get("mutation_rate", DEFAULT_MUTATION_RATE)
    recombination_rate = parameter.get("recombination_rate", DEFAULT_RECOMBINATION_RATE)
    
    #Run the model TODO:change silent to true for actual runs
    Main.main(num_clusters=numClusters, migration_rates_modifier=m, population_modifier=pop, 
              mutation_rate=mutation_rate, recombination_rate=recombination_rate, silent=True)
    
    
    #Read in the output data
    outDict = {}
    
    for year in ["2015", "2019", "2023"]:
        with open(Path(f"../data/Output_Data/diversities_{year}.csv"), mode='r', newline='', encoding='utf-8') as csvfile:
            div2015 = csv.DictReader(csvfile)
            diversities_list = []
            for row in div2015:
                value = next(iter(row.values()))
                if value is not None and value.strip() != "":
                    diversities_list.append(float(value.strip()))
            diversities = np.array(diversities_list)
            outDict[f"{year}_diversity"] = diversities
        
        with open(Path(f"../data/Output_Data/divergences_{year}.csv"), mode='r', newline='', encoding='utf-8') as csvfile:
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
        "numClusters": prior_distributions["numClusters"].rvs(),
        "mutation_rate": prior_distributions["mutation_rate"].rvs(),
        "recombination_rate": prior_distributions["recombination_rate"].rvs()
    }


def run_abc_simulation_sensitivity(output_csv="../out/abc_results.csv"):
    '''
    Run ABC simulations with sensitivity analysis: 35 iterations in 5 sets of 7.
    Each set varies a single parameter according to its distribution while keeping others fixed.
    
    Parameter sets:
    1. Vary migration rate (m)
    2. Vary population size (pop)
    3. Vary number of clusters (numClusters)
    4. Vary mutation rate (with log-normal distribution)
    5. Vary recombination rate (with log-normal distribution)
    
    :param output_csv: Path to the output CSV file
    '''
    
    observed_data = getObservedData()
    
    # Determine if we need to write the header
    csv_exists = Path(output_csv).exists()
    
    # Define CSV columns
    fieldnames = ["iteration", "set", "varied_parameter", "m", "pop", "numClusters", 
                  "mutation_rate", "recombination_rate",
                  "diversity_loss_2015", "diversity_loss_2019", "diversity_loss_2023",
                  "divergence_loss_2015", "divergence_loss_2019", "divergence_loss_2023",
                  "total_loss"]
    
    # Get baseline parameters (sample once to establish central values for fixed parameters)
    baseline = {
        "m": prior_distributions["m"].rvs(),
        "pop": prior_distributions["pop"].rvs(),
        "numClusters": prior_distributions["numClusters"].rvs(),
        "mutation_rate": DEFAULT_MUTATION_RATE,
        "recombination_rate": DEFAULT_RECOMBINATION_RATE
    }
    
    # Define the 5 parameter sets
    parameter_sets = [
        {"name": "migration_rate", "param_key": "m", "values": [prior_distributions["m"].rvs() for _ in range(7)]},
        {"name": "population_size", "param_key": "pop", "values": [prior_distributions["pop"].rvs() for _ in range(7)]},
        {"name": "num_clusters", "param_key": "numClusters", "values": [prior_distributions["numClusters"].rvs() for _ in range(7)]},
        {"name": "mutation_rate", "param_key": "mutation_rate", "values": [prior_distributions["mutation_rate"].rvs() for _ in range(7)]},
        {"name": "recombination_rate", "param_key": "recombination_rate", "values": [prior_distributions["recombination_rate"].rvs() for _ in range(7)]}
    ]
    
    iteration = 0
    
    with open(output_csv, mode='a', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        # Write header only if file is new
        if not csv_exists:
            writer.writeheader()
        
        # Iterate through each parameter set
        for set_idx, param_set in enumerate(parameter_sets):
            set_name = param_set["name"]
            param_key = param_set["param_key"]
            values = param_set["values"]
            
            print(f"\n=== Set {set_idx + 1}: Varying {set_name} ===")
            
            for i, value in enumerate(values):
                # Create parameters for this iteration
                parameters = baseline.copy()
                parameters[param_key] = value
                
                print(
                    f"Running iteration {iteration + 1}/35 (Set {set_idx + 1}, Sample {i + 1}/7) "
                    f"varying {set_name}: {param_key}={value:.6g}..."
                )
                
                try:
                    # Run the model
                    simulated_data = model(parameters)
                    
                    # Calculate losses
                    losses = calculate_losses(observed_data, simulated_data)
                    
                    # Prepare row for CSV
                    row = {
                        "iteration": iteration,
                        "set": set_idx + 1,
                        "varied_parameter": set_name,
                        "m": parameters["m"],
                        "pop": parameters["pop"],
                        "numClusters": parameters["numClusters"],
                        "mutation_rate": parameters["mutation_rate"],
                        "recombination_rate": parameters["recombination_rate"],
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
                    iteration += 1
                    
                except Exception as e:
                    print(f"  Error in iteration {iteration}: {e}")
                    iteration += 1
                    continue
    
    print(f"\n=== Sensitivity Analysis Complete ===")
    print(f"Total iterations: {iteration}/35. Results saved to {output_csv}")


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
    fieldnames = ["iteration", "m", "pop", "numClusters", "mutation_rate", "recombination_rate",
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
                f"numClusters={parameters['numClusters'] * 33}, "
                f"mutation_rate={parameters['mutation_rate']:.6g}, "
                f"recombination_rate={parameters['recombination_rate']:.6g}..."
            )
            
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
                    "mutation_rate": parameters["mutation_rate"],
                    "recombination_rate": parameters["recombination_rate"],
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
    # Run 35 iterations of ABC sampling with sensitivity analysis
    # Each of 5 sets varies a single parameter (7 samples per set)
    run_abc_simulation_sensitivity(output_csv="../out/abc_results.csv")