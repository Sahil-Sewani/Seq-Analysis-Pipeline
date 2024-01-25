"""
This program functions to make individual fasta files for each 1:1 ortholog pair between all species pairs.
It does this by creating directories if they do not exist and skipping them if they already exist.
It processes amino acid files first and then the DNA file.

The 1:1 ortholog file needs to be written into the code at the top of the pair_species_genes method.
"""

import os

# Method to read the given species pair file and put them into a list of tuples
def pairs_of_eukaryotic_species(user_input_data='oma-pair_eukaryote.txt'):
    speciesPairs = []
    with open(user_input_data, 'r') as user_input_species_list:
        for line in user_input_species_list:
            if line.strip() != "":
                new_line = line.split()
                speciesPairs.append((new_line[0], new_line[1]))
    return speciesPairs

# Method to find 1:1 orthologs for given species pairs and organize them into a dictionary
def ortholog_finder(list_of_species_pairs):
    filtered = 'oma-pair_eukaryote.txt'
    pairs_and_lines = {pair: [] for pair in list_of_species_pairs}

    with open(filtered, 'r') as orthopairs:
        counter = 0
        for line in orthopairs:
            counter += 1
            if counter % 10000000 == 0:
                print("Processing line {:,} of ortholog pairs".format(counter))
            new_line = line.split()
            sp1 = new_line[0][0:5]
            sp2 = new_line[1][0:5]

            if (sp1, sp2) in list_of_species_pairs:
                pairs_and_lines[(sp1, sp2)].append((new_line[0][5:].strip(), new_line[1][5:].strip()))
            elif (sp2, sp1) in list_of_species_pairs:
                pairs_and_lines[(sp2, sp1)].append((new_line[1][5:].strip(), new_line[0][5:].strip()))

    return pairs_and_lines

# Method to create individual fasta files for each species pair
def pair_species_genes(ortholog_database):
    for key, item in ortholog_database.items():
        print("Writing output for species pair {} / {}".format(key[0], key[1]))

        species1_file = 'species_fasta/' + key[0] + '_aa' + '.fasta'
        species2_file = 'species_fasta/' + key[1] + '_aa' + '.fasta'

        folderDirectory = ('pairs_fasta/' + key[0] + '_' + key[1])
        if not os.path.isdir(folderDirectory):
            os.makedirs(folderDirectory)

        FirstSpecies = {}
        SecondSpecies = {}
        with open(species1_file, 'r') as firstFile:
            for line in firstFile:
                sp1key = line[7:].strip()
                sequence = firstFile.readline()
                FirstSpecies[sp1key] = sequence

        with open(species2_file, 'r') as secondFile:
            for line in secondFile:
                sp2key = line[7:].strip()
                sequence = secondFile.readline()
                SecondSpecies[sp2key] = sequence

        with open(folderDirectory + '/' + key[0] + '_' + key[1] + '_pairs.txt', "w") as outfile:
            for x in item:
                outfile.write(key[0] + x[0] + "\t" + key[1] + x[1] + "\n")

        counter = 1
        for x in item:
            padNum = str(counter).rjust(5, '0')
            outputFileName = folderDirectory + '/' + key[0] + '_' + key[1] + '_aa' + '_' + padNum + '.fasta'

            with open(outputFileName, 'w') as output:
                output.write('> ' + key[0] + x[0] + '\n')
                output.write(FirstSpecies[x[0]] + '\n')
                output.write('> ' + key[1] + x[1] + '\n')
                output.write(SecondSpecies[x[1]] + '\n')
            counter += 1

# Main function to orchestrate the whole process
def main():
    # Step 1: Get the name of the file containing the species pairs
    userdata = input("Enter the file of species pairs: ")

    # Step 2: Read in the file
    species_pairs_list = pairs_of_eukaryotic_species(userdata)

    # Step 3: Call a function to parse the ortholog pairs file and find matching orthologs for the pairs fetched in Step 2
    ortholog_database = ortholog_finder(species_pairs_list)

    # Step 4: For every species pair found in Step 2, create the folder and all the files for that pair
    pair_species_genes(ortholog_database)

if __name__ == "__main__":
    main()

