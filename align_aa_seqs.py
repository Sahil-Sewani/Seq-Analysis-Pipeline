"""
This program aligns a collection of paired amino acid sequence files using the MUSCLE alignment tool.
"""

import glob
import subprocess
from tqdm import tqdm  # tqdm is a library for adding progress bars to loops
from datetime import datetime

# Method to read the given species pair file and put them into a list of tuples
def pairs_of_eukaryotic_species(user_input_data='oma-pair_eukaryote.txt'):
    speciesPairs = []
    with open(user_input_data, 'r') as user_input_species_list:
        for line in user_input_species_list:
            new_line = line.split()
            speciesPairs.append((new_line[0], new_line[1]))
    return speciesPairs

# Method to align pairs of amino acid sequence files for a given species pair
def align_pairs(pair):
    sp1 = pair[0]
    sp2 = pair[1]
    file_pattern = 'pairs_fasta/' + sp1 + '_' + sp2 + "/*_aa_?????.fasta"
    filelist = glob.glob(file_pattern)
    print("Working on {} / {}".format(sp1, sp2))
    for f in tqdm(filelist):  # tqdm is an automatic progress indicator
        outf = f[:f.find(".fasta")] + "_aligned.fasta"
        command = ["muscle", "-align", f, "-output", outf]
        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
    print()

# Main function to orchestrate the whole process
def main():
    start = datetime.now()
    # Step 1: Get the name of the file containing the species pairs
    userdata = input("Enter the file of species pairs: ")

    # Step 2: Read in the file
    species_pairs_list = pairs_of_eukaryotic_species(userdata)

    # Step 3: Align all files for each pair from user data
    for pair in species_pairs_list:
        align_pairs(pair)

    print("Total run time:", datetime.now() - start)

if __name__ == "__main__":
    main()
