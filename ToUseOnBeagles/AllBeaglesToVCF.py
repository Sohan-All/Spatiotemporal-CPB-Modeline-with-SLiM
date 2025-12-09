import subprocess
import os

def convertAll():
    for i in range(1, 18): 
        beagleZipped = f"cpbWGS_genolike_chr{i}.cpbWGS_genolike_chr{i}.beagle.gz.phased.gz"
        beagleUnZipped = f"cpbWGS_genolike_chr{i}.beagle"
        outputVCF = f"chr{i}_cpb.vcf"

        if os.path.exists(outputVCF + ".gz"):
            print(f"chromosome {i} VCF already exists")
            continue

        print(f"Started converting chromosome {i} to VCF. ")
        print(f"Unzipping beagle for chromosome {i}. ")
        subprocess.run(['gunzip', '-c', beagleZipped], stdout=open(beagleUnZipped, 'wb'))
        print("Converting chromosome {i} beagle file to VCF. ")
        subprocess.run(['python', 'ConvertBeagleToVCF.py', beagleUnZipped, outputVCF])
        subprocess.run(['rm', beagleUnZipped])
        print(f"Zipping VCF for chromosome {i}.")
        subprocess.run(['bgzip', outputVCF])
        #subprocess.run(['rm', outputVCF])
        print(f"Finished for chromosome {i}.")
    
if __name__ == "__main__":
    convertAll()
