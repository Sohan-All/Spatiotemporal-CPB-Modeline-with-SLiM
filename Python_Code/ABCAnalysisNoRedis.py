import csv
import numpy as np
import sys
import shutil

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


def read_parameters_from_csv(csv_path):
    '''
    Read parameter configurations from a CSV file.
    
    Expected CSV columns: m, pop, numClusters, mutation_rate, recombination_rate
    Each row represents one simulation to run.
    
    :param csv_path: Path to the CSV file with parameters
    :return: List of dictionaries, each containing parameters for one simulation
    '''
    parameters_list = []
    
    try:
        with open(csv_path, mode='r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            
            if reader.fieldnames is None:
                raise ValueError(f"CSV file {csv_path} is empty or has no headers")
            
            # Validate that all required columns are present
            required_cols = {"m", "pop", "numClusters", "mutation_rate", "recombination_rate"}
            csv_cols = set(reader.fieldnames)
            missing_cols = required_cols - csv_cols
            
            if missing_cols:
                raise ValueError(f"CSV file missing required columns: {missing_cols}. "
                                f"Required columns: {required_cols}")
            
            for row_idx, row in enumerate(reader, start=2):  # start=2 because row 1 is header
                try:
                    parameters = {
                        "m": float(row["m"]),
                        "pop": int(float(row["pop"])),  # Convert to float first to handle scientific notation
                        "numClusters": int(float(row["numClusters"])),
                        "mutation_rate": float(row["mutation_rate"]),
                        "recombination_rate": float(row["recombination_rate"])
                    }
                    parameters_list.append(parameters)
                except ValueError as e:
                    print(f"Warning: Row {row_idx} in {csv_path} has invalid values: {e}")
                    continue
        
        if not parameters_list:
            raise ValueError(f"No valid parameter configurations found in {csv_path}")
        
        print(f"Loaded {len(parameters_list)} parameter configuration(s) from {csv_path}")
        return parameters_list
    
    except FileNotFoundError:
        raise FileNotFoundError(f"Input CSV file not found: {csv_path}")


def run_sims_from_csv(input_csv, output_csv="../out/abc_results.csv", simToRun=-1):
    '''
    Run ABC simulations with parameters specified in a CSV file.
    Each row in the input CSV represents one simulation.
    Detailed simulation outputs (diversities and divergences) are saved to detailed_sim_results folder.
    
    :param input_csv: Path to the CSV file with input parameters
    :param output_csv: Path to the output CSV file for results
    :param simToRun: Index of the specific simulation to run (if -1, run all)
    '''
    
    try:
        parameters_list = read_parameters_from_csv(input_csv)
    except Exception as e:
        print(f"Error reading input CSV: {e}")
        return
    
    observed_data = getObservedData()
    
    # Determine if we need to write the header
    csv_exists = Path(output_csv).exists()
    
    # Create detailed results directory
    output_dir = Path(output_csv).parent
    detailed_results_dir = output_dir / "detailed_sim_results"
    detailed_results_dir.mkdir(parents=True, exist_ok=True)
    print(f"Detailed results will be saved to: {detailed_results_dir}")
    
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
        
        for iteration, parameters in enumerate(parameters_list):
            if simToRun != -1 and iteration != simToRun:
                continue

            print(
                f"Running iteration {iteration + 1}/{len(parameters_list)} "
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
                
                # Copy detailed results for this iteration
                iteration_dir = detailed_results_dir / f"run{iteration + 1}"
                iteration_dir.mkdir(parents=True, exist_ok=True)
                
                for year in ["2015", "2019", "2023"]:
                    # Copy diversities
                    diversity_src = Path(f"../data/Output_Data/diversities_{year}.csv")
                    diversity_dest = iteration_dir / f"diversities_{year}.csv"
                    if diversity_src.exists():
                        shutil.copy2(diversity_src, diversity_dest)
                    
                    # Copy divergences
                    divergence_src = Path(f"../data/Output_Data/divergences_{year}.csv")
                    divergence_dest = iteration_dir / f"divergences_{year}.csv"
                    if divergence_src.exists():
                        shutil.copy2(divergence_src, divergence_dest)
                
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
                print(f"  Detailed results saved to: {iteration_dir}")
                
            except Exception as e:
                print(f"  Error in iteration {iteration}: {e}")
                continue
    
    print(f"\n=== CSV-based simulations complete ===")
    print(f"Total iterations: {len(parameters_list)}. Results saved to {output_csv}")
    print(f"Detailed simulation data saved to {detailed_results_dir}")



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
    if len(sys.argv) != 2:
        print("Usage: python ABCAnalysisNoRedis.py <integer>")
        sys.exit(1)
    
    simToRun = int(sys.argv[1])
    run_sims_from_csv("./sample_inputs.csv", "../out/abc_results.csv", simToRun)
    