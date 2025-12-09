import sys


def main():    
    if len(sys.argv) != 3:
        print("Usage: python ConvertBeagleToVCF.py <input_specifierMatrix_file> <output_populations_file>")
        sys.exit(1)
        
    input_specifierMatrix_file = sys.argv[1]
    output_populations_file = sys.argv[2]
    with open(input_specifierMatrix_file, 'r') as fin, open(output_populations_file, 'w') as fout:
        for line in fin:
            parts = line.strip().split(",")
            populationName = parts[0]
            for i in range(3, len(parts)):
                if parts[i] == "":
                    continue
                id = "S" + str(int(float(parts[i])))
                fout.write(f"{id}\t{populationName}\n")
                
    
main()