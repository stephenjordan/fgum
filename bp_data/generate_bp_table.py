import re
import math

# Initialize four separate lists
k_list = []
D_list = []
threshold_list = []
uncertainty_list = []
threshmax_list = []
threshmin_list = []

# Regular expressions to catch the relevant lines in the text file
# Matches lines like: Executing with k=3, D=4...
kd_pattern = re.compile(r"Executing with k=(\d+),\s*D=(\d+)\.\.\.")
# Matches lines like: threshold appears to be at least 0.167055845293798 - 0.167055964503087
bounds_pattern = re.compile(r"threshold appears to be at least ([\d\.]+) - ([\d\.]+)")

# Open and parse the file
with open('de_results.txt', 'r') as file:
    current_k = None
    current_D = None
    
    for line in file:
        # Check if the line contains k and D
        kd_match = kd_pattern.search(line)
        if kd_match:
            current_k = int(kd_match.group(1))
            current_D = int(kd_match.group(2))
        
        # Check if the line contains the threshold bounds
        bounds_match = bounds_pattern.search(line)
        if bounds_match and current_k is not None and current_D is not None:
            lower_bound = float(bounds_match.group(1))
            upper_bound = float(bounds_match.group(2))
            
            # Compute threshold and uncertainty
            threshold = (lower_bound + upper_bound) / 2.0
            uncertainty = abs(upper_bound - lower_bound)
            
            # Append to the lists
            k_list.append(current_k)
            D_list.append(current_D)
            threshold_list.append(threshold)
            uncertainty_list.append(uncertainty)
            threshmax_list.append(upper_bound)
            threshmin_list.append(lower_bound)
            
            # Reset current k and D to ensure we don't duplicate data if the file formatting gets weird
            current_k = None
            current_D = None

max_uncertainty = threshmax_list[0] - threshmin_list[0]

for i in range(1,len(threshmax_list)):
    uncertainty = threshmax_list[i] - threshmin_list[i]
    if uncertainty > max_uncertainty:
        max_uncertainty = uncertainty

print("max uncertainty in BP threshold: " + str(max_uncertainty))

#semicircle law
def s(thresh):
    return 0.5+math.sqrt(thresh*(1.0-thresh))

print("k", "D", "s", sep="\t")
for i in range(len(k_list)):
    print(f"{k_list[i]}\t{D_list[i]}\t{s(threshmin_list[i]):.4f}")