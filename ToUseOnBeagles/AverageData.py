import subprocess
import os
import sys

def main():    
    numChrs = 17
    
    for year in {"2015", "2019", "2023"}:
        #Get indices of each site to ensure the sites are in the same order for each statistic
        siteNames = []
        with open(f"specifier_matrix_{year}.csv", 'r') as f:
            lines = f.readlines()
            for line in lines:
                values = line.strip().split(',')
                siteNames.append(values[0])
            
        numSites = len(siteNames)
                    
        piData = [[] for _ in range(numSites)] 
        dxyData = [[[] for _ in range(numSites)] for _ in range(numSites)] 
        fstData = [[[] for _ in range(numSites)] for _ in range(numSites)] 
        
        for i in range(1, numChrs + 1): 
            folderName = f"statsChr{i}_{year}"
            piPath = os.path.join(folderName, f"chr{i}_{year}_pi.txt")
            dxyPath = os.path.join(folderName, f"chr{i}_{year}_dxy.txt")
            fstPath = os.path.join(folderName, f"chr{i}_{year}_fst.txt")
            
            
            with open(piPath, 'r') as f:
                data = f.readlines()
                for line in data:
                    if line.startswith("pop"):
                        continue
                    values = line.strip().split()
                    siteIdx = siteNames.index(values[0])
                    piData[siteIdx].append(float(values[4]))
                    
            with open(dxyPath, 'r') as f:
                data = f.readlines()
                for line in data:
                    if line.startswith("pop"):
                        continue
                    values = line.strip().split()
                    siteIdx1 = siteNames.index(values[0])
                    siteIdx2 = siteNames.index(values[1])
                    dxyData[siteIdx1][siteIdx2].append(float(values[5]))
                    dxyData[siteIdx2][siteIdx1].append(float(values[5]))
            
            with open(fstPath, 'r') as f:
                data = f.readlines()
                for line in data:
                    if line.startswith("pop"):
                        continue
                    values = line.strip().split()
                    siteIdx1 = siteNames.index(values[0])
                    siteIdx2 = siteNames.index(values[1])
                    fstData[siteIdx1][siteIdx2].append(float(values[5]))
                    fstData[siteIdx2][siteIdx1].append(float(values[5]))
            
        
        averagedPis = [sum(piData[i])/len(piData[i]) for i in range(numSites)]
        averagedDxy = [[0 if k == i else sum(dxyData[i][k])/len(dxyData[i][k]) for k in range(numSites)] for i in range(numSites)]
        averagedFst = [[0 if k == i else sum(fstData[i][k])/len(fstData[i][k]) for k in range(numSites)] for i in range(numSites)]
        
        with open(f"finalStats/averaged_pi_{year}.csv", 'w') as f:
            for site, pi in zip(siteNames, averagedPis):
                f.write(f"{pi}\n")
                
        with open(f"finalStats/averaged_dxy_{year}.csv", 'w') as f:
            for i in range(numSites):
                line = ""
                for k in range(numSites):
                    line += f"{averagedDxy[i][k]},"
                f.write(line[:-1] + "\n")
        
        with open(f"finalStats/averaged_fst_{year}.csv", 'w') as f:
            for i in range(numSites):
                line = ""
                for k in range(numSites):
                    line += f"{averagedFst[i][k]},"
                f.write(line[:-1] + "\n")
                
                
if __name__ == "__main__":
    main()