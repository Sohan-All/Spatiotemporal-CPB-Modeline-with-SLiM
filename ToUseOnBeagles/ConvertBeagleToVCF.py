import datetime

def generate_VCF_from_beagle(beagleFilepath="../../../Data for modeling/cpbWGS_genolike_chr9.cpbWGS_genolike_chr9.beagle.gz.phased", outputVCFpath="../data/Genetic_Data/chr9CPB.vcf"):
    """
    Generate a vcf file from beagle genetic data.
    This function reads the genetic data and creates a VCF file.
    
    Parameters:
    beagleFilepath (str): Path to the beagle genetic data file. This file is not in the github repository and should be provided separately.
    outputVCFpath (str): Path where the output VCF file will be saved.
    
    Returns:
    None
    """
    CUTOFF = 10000
    base_map = {'0': 'A', '1': 'C', '2': 'G', '3': 'T'}
    
    
    print("started creating vcf")
    
    with open(beagleFilepath, 'r') as fin:
        with open(outputVCFpath, "w") as fout:
        
            # skip the first/header line
            next(fin, None)
            # write the VCF header
            num_samples = (len(fin.readline().strip().split()) - 2) // 2
            fout.write(generate_VCF_header(num_samples))
            
            counter = 0
            for line in fin:
                counter += 1
                if counter > CUTOFF:
                    break
                parts = line.strip().split()
                
                CHROM = 9
                POS = int(parts[1].split('_')[1])
                ID = "."
                REF = base_map.get(parts[2], 'N') #reference base is the first allele's base at this position
                QUAL = "100"
                FILTER = "PASS"
                INFO = "."
                FORMAT = "GT"
                
                alts = []
                for i in range(2, len(parts)):
                    if parts[i] != parts[2]:
                        alt_base = base_map.get(parts[i])
                        if alt_base not in alts:
                            alts.append(alt_base)
                            
                ALT = ",".join(alts)
                
                samples = []
                
                for i in range(2, len(parts), 2):
                    allele1 = base_map.get(parts[i])
                    allele2 = base_map.get(parts[i+1])
                    strToAdd = ""
                    if allele1 == REF:
                        strToAdd += "0"
                    else:
                        strToAdd += f"{alts.index(allele1) + 1}"
                    strToAdd += "|"
                    if allele2 == REF:
                        strToAdd += "0"
                    else:
                        strToAdd += f"{alts.index(allele2) + 1}"
                        
                    samples.append(strToAdd)
                    
                fout.write(f"{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}\t" + "\t".join(samples) + "\n")
                        
    print("finished creating vcf")
    
    
def generate_VCF_header(num_samples):
    """
    Generate the header for a VCF file.
    
    Parameters:
    num_samples (int): The number of samples in the VCF file.
    
    Returns:
    str: The VCF header as a string.
    """
    header = "##fileformat=VCFv4.5\n"
    header += f"##datecreated={datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")}\n"
    header += "##source=ConvertBeagleToFasta.py\n"
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    
    for i in range(1, num_samples + 1):
        header += f"\tS{i}"
    
    header += "\n"
    return header
            
