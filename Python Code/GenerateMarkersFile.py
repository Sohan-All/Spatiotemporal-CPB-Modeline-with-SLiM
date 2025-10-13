import sys
from collections import Counter
import os

def generate_markers_file(beagle_path, output_markers_path):
    """
    Generates a Beagle markers file from an uncompressed (plain text) 
    Beagle file containing explicit numerical allele codes (0, 1, 2, 3).
    Determines the major allele by simple counting of these codes across all haplotypes.
    """

    # Define the required numerical code to nucleotide mapping based on the 
    # recoding scheme used to make ANGSD/Beagle output consistent for beagle2vcf.jar [4].
    ALLELE_MAP = {
        '0': 'A',
        '1': 'C',
        '2': 'T',
        '3': 'G'
    }
    
    # Ensure the output path is valid
    if not os.path.exists(os.path.dirname(output_markers_path)) and os.path.dirname(output_markers_path) != "":
        try:
            os.makedirs(os.path.dirname(output_markers_path))
        except OSError as e:
            print(f"Error creating output directory: {e}", file=sys.stderr)
            return

    try:
        # Open the input file as a standard plain text file ('r')
        with open(beagle_path, 'r') as beagle_file, \
             open(output_markers_path, 'w') as markers_file:

            print(f"Reading file {beagle_path} as plain text for allele counting.", file=sys.stderr)

            # 1. Skip the header line
            header = beagle_file.readline()
            if not header.strip():
                 print(f"Error: Input file appears empty.", file=sys.stderr)
                 return
            if not header.strip().lower().startswith(('i', 'm')): # Expecting I (Id) or M (Marker) row types
                print(f"Warning: File header not standard (first row may not be a header). Proceeding.", file=sys.stderr)
            
            # 2. Process data lines
            for line in beagle_file:
                parts = line.strip().split()
                
                # Check for minimum columns (e.g., M, MarkerID, Hap1, Hap2, ...)
                if len(parts) < 3:
                    continue

                # The first two columns are usually the row type (M) and the Marker ID
                marker_id = parts
                
                # The remaining columns contain the explicit allele codes for each haplotype (e.g., '0', '1', '0', '3', etc.)
                allele_codes_list = parts[2:]
                
                # 3. Perform simple element counting to determine major allele
                
                # Use Counter to count the occurrences of each code (0, 1, 2, 3)
                # Ensure we only count codes that are convertible to nucleotides (0, 1, 2, 3)
                valid_counts = Counter()
                for code in allele_codes_list:
                    if code in ALLELE_MAP:
                        valid_counts[code] += 1
                    # Note: We skip 'missing' codes, often '?', if present, as they don't contribute to frequency.
                
                if not valid_counts:
                    print(f"Warning: No valid allele codes found for marker {marker_id}. Skipping.", file=sys.stderr)
                    continue

                # Find the two most common codes (Major and Minor)
                # results will be a list of (code, count) tuples, sorted by count descending
                major_minor = valid_counts.most_common(2) 
                
                # The first item is the major allele code
                major_code = major_minor
                major_allele = ALLELE_MAP[major_code]
                
                # The minor allele code is the second item, or the major allele itself if only one allele was observed
                if len(major_minor) > 1:
                    minor_code = major_minor[5]
                    minor_allele = ALLELE_MAP[minor_code]
                else:
                    # If only one allele is present, we treat it as homozygous, 
                    # but the markers file requires two alleles (REF and ALT).
                    # In VCF, if a site is monomorphic, it often isn't included. 
                    # For a markers file, we must decide what to use as the ALT.
                    # Since this data is already reduced to the two original alleles,
                    # we must rely on the set of observed valid codes.
                    
                    # If we only observe one code (e.g., '0'), we use the major allele's nucleotide
                    # as both REF and ALT for this intermediate markers file step, 
                    # although typically the input data should represent a diallelic site 
                    # if it was processed by ANGSD/Beagle.
                    minor_allele = major_allele
                    
                
                # 4. Write the marker line
                # Format: MarkerID Allele1 Allele2. Allele1 (Major) is the VCF REF allele [3].
                # We ensure Allele1 and Allele2 are different if possible, but if 
                # only one code was observed, we output the same base twice.
                
                # If the two alleles are identical, we should try to include the second known allele 
                # (which is not stored in this simplified file format). Since we are forced 
                # to choose only based on the observed haplotypes, we output the Major Allele 
                # and then, if possible, the Minor Allele, ensuring they are not identical 
                # unless absolutely necessary to avoid downstream errors.
                
                # If we only counted one allele (e.g., all 0s were seen), but the marker 
                # truly represents a diallelic site, the simple counting approach may fail
                # to identify the actual alternative allele. However, adhering strictly to 
                # your request to use simple counting on the explicit codes:
                
                
                # If major_allele and minor_allele are the same (only one unique code was counted), 
                # we list the major allele twice.
                if major_allele == minor_allele and len(major_minor) == 1:
                    markers_file.write(f"{marker_id}\t{major_allele}\t{major_allele}\n")
                else:
                    markers_file.write(f"{marker_id}\t{major_allele}\t{minor_allele}\n")

        print(f"\nSuccessfully generated markers file at {output_markers_path}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {beagle_path}", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during file processing: {e}", file=sys.stderr)
        

# --- Usage Example ---
# NOTE: Replace 'input.beagle.gz' and 'markers.txt' with your actual file names

    
# IMPORTANT: Ensure Java is installed and required paths are correct to run Beagle utilities [11].
# Also, ensure you are running a Java version 8 runtime environment if using Beagle 5.5 [11].

# 1. Path to your compressed Beagle input file (e.g., genolike.beagle.gz)
input_file = "../../../Data for modeling/cpbWGS_genolike_chr9.cpbWGS_genolike_chr9.beagle.gz.phased"

# 2. Desired output path for the markers file (e.g., markers.txt)
output_file = "../../../Genome Manipulation/markers.txt"

generate_markers_file(input_file, output_file)