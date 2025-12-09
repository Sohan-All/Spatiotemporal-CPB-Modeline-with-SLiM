import subprocess
import os
import sys

def main():
    for year in {"2015", "2019", "2023"}:
        for i in range(1, 18): 
            inVCF = f"chr{i}_cpb.vcf.gz"
            popFile = "popFile" + year
            outFolder = f"statsChr{i}_{year}"
            prefix = f"chr{i}_{year}"
            print(f"calculating pi for chr {i} year {year}") 
            os.makedirs(outFolder, exist_ok=True)
            
            with open("temp.txt", 'w') as output_file:
                cmd = f"zcat {inVCF} | tail -n 1"
                subprocess.run(cmd, shell=True, stdout=output_file, check=True)

            ending_base = 0
            # read temp.txt's second line and take the second element
            with open('temp.txt', 'r') as f:
                lines = f.readlines()
            
                line = lines[0].strip()
                parts = line.split()
                ending_base = parts[1]
            
            subprocess.run(['pixy', '--stats', 'pi', 'dxy', 'fst', '--vcf', inVCF, '--populations', popFile, "--window_size", ending_base, "--n_cores", "8", "--output_folder", outFolder, "--output_prefix", prefix, "--bypass_invariant_check", "--silent"])
    
if __name__ == "__main__":
    main()
