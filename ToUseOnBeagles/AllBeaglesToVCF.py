import subprocess

def convertAll():
    for i in range(1, 18):
        beagleZipped = f"cpbWGS_genolike_chr{i}.cpbWGS_genolike_chr{i}.beagle.gz.phased.gz"
        beagleUnZipped = f"cpbWGS_genolike_chr{i}.beagle"
        outputVCF = f"chr{i}_cpb.vcf"
        subprocess.run(['gunzip', '-c', beagleZipped], stdout=open(beagleUnZipped, 'wb'))
        subprocess.run(['python', 'ConvertBeagleToVCF.py', beagleUnZipped, outputVCF])
        subprocess.run(['rm', beagleUnZipped])
        subprocess.run(['bgzip', outputVCF])
        subprocess.run(['rm', outputVCF])
        print(f"Finished converting chromosome {i} to VCF.")
    
if __name__ == "__main__":
    convertAll()