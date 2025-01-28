import glob
import os
import re
import sys

## How to run this program?
# python Summarize_top_targets.py safety90 AllDatasets_EachPatient500Cells_optimal_doublets
# python Summarize_top_targets.py safety90 AllDatasets_EachPatient500Cells_optimal_triplets
# python Summarize_top_targets.py safety90 AllDatasets_EachPatient500Cells_optimal_quadruples
# python Summarize_top_targets.py safety90 AllDatasets_EachPatient500Cells_optimal_quintuplets

# Function to remove prefix and suffix
def extract_substring(filename, prefix, suffix_pattern):
    if filename.startswith(prefix):
        trimmed = filename[len(prefix):]
    trimmed = re.sub(suffix_pattern, "", trimmed)
    return trimmed


def remove_duplicates(data):
    # Create a dictionary to hold unique elements based on the third and fourth elements
    unique_elements = {}
    for item in data:
        # Convert the third element (list of strings) to a tuple to make it hashable
        key = (tuple(item[2]), item[3])
        if key not in unique_elements:
            unique_elements[key] = item
    # Return the values of the dictionary as a list
    return list(unique_elements.values())


if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python script.py <folder> <prefix>")
        print("Example 1: python script.py Stefan_AA_11_allTumorCells BRCA_GSE176078_optimal_triplets")
        print("Example 2: python script.py '' BRCA_GSE176078_optimal_triplets (for ./output)")
        sys.exit(1)

    # Define the folder and prefix
    folder = sys.argv[1]  # e.g., Stefan_AA_11_allTumorCells or empty string
    prefix = sys.argv[2]  # e.g., BRCA_GSE176078_optimal_triplets

    # Construct the directory path
    if folder:  # If folder is not an empty string
        directory = os.path.join('./output', folder)
    else:  # If folder is an empty string, use ./output directly
        directory = './output'

    # Check if the directory exists
    if not os.path.exists(directory):
        raise FileNotFoundError(f"The directory '{directory}' does not exist.")

    # Define the suffix pattern for matching files
    suffix_pattern = r"_\d+\.txt"

    # Construct the output file path
    fnOut = os.path.join(directory, f"{prefix}_merged.txt")

    # Print confirmation
    print(f"Processing files in directory: {directory}")
    print(f"Merged output will be saved to: {fnOut}")

    # Use glob to get all files with the specified prefix
    file_pattern = os.path.join(directory, f'{prefix}*')
    files = glob.glob(file_pattern)
    files = [file for file in files if 'GAlog' not in os.path.basename(file)]
    files = [file for file in files if re.search(suffix_pattern, os.path.basename(file))]

    # Loop through each file and read it
    all_CAR_combinations = []

    for file in files:
        basename = os.path.basename(file)
        gate_name = extract_substring(basename, prefix, suffix_pattern)[1:]
        CAR_combinations = open(file, 'r').readlines()[1:]
        for i in range(len(CAR_combinations)):  # different predicted gene combinations
            gene_temp_str,score_str = CAR_combinations[i].split('\t')[1:]
            score = float(score_str)
            gene_temp = gene_temp_str[2:-2].split("', '")
            gene_temp = sorted(gene_temp)
            all_CAR_combinations.append([gate_name, gene_temp_str, gene_temp, score])
    # Remove duplicates
    all_CAR_combinations = remove_duplicates(all_CAR_combinations)
    all_CAR_combinations.sort(key=lambda x: -x[3])
    for c in all_CAR_combinations[:50]:
        print(c[0:2], round(c[3], 3))

    fhOut = open(fnOut, "w")
    fhOut.write("LogicGates\tGenes\tScore\n")
    for c in all_CAR_combinations:
        content = "\t".join([c[0], c[1], str(c[3])])
        fhOut.write(content + "\n")
    fhOut.close()
