import re
from collections import defaultdict

def parse_max_s_data(file_path):
    # Dictionary to store the maximum s found for each (k, D) pair
    # Key: (k, D), Value: {'max_s': int, 'm': int, 'sweeps': str}
    max_results = defaultdict(lambda: {'max_s': 0, 'm': 0, 'sweeps': "1E4"})
    
    content_re = re.compile(r"(\d+) clauses satisfied out of (\d+)")
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or ':' not in line:
                    continue
                
                filename, content = line.split(':', 1)
                
                # 1. Infer k and D from the filename logic
                file_parts = filename.replace('.txt', '').split('_')
                k = int(file_parts[-4]) # Adjusted for the trailing trial index (_1, _2, etc)
                D = int(file_parts[-3])
                
                # 2. Extract s and m
                match = content_re.search(content)
                if match:
                    s = int(match.group(1))
                    m = int(match.group(2))
                    
                    # 3. Update the entry only if this trial's s is the new maximum
                    if s > max_results[(k, D)]['max_s']:
                        max_results[(k, D)]['max_s'] = s
                        max_results[(k, D)]['m'] = m
                        
        # Convert dictionary to a sorted list for display
        sorted_data = []
        for (k, D), val in max_results.items():
            sorted_data.append({
                'k': k,
                'D': D,
                's': val['max_s'],
                'm': val['m'],
                'sweeps': "1E6" # Overriding to 1E6 as per your instruction
            })
        
        sorted_data.sort(key=lambda x: (x['k'], x['D']))
        return sorted_data

    except FileNotFoundError:
        print(f"Error: {file_path} not found.")
        return []

# Run the parser
data = parse_max_s_data('summary.txt')

# Display the grouped results
print(f"{'k':<4} {'D':<4} {'Sweeps':<8} {'Max s':<8} {'m':<8} {'Best Ratio':<12}")
print("-" * 55)
for entry in data:
    ratio = f"{entry['s']/entry['m']:.5f}"
    print(f"{entry['k']:<4} {entry['D']:<4} {entry['sweeps']:<8} {entry['s']:<8} {entry['m']:<8} {ratio:<12}")