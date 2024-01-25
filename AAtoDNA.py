import os
import glob
import tqdm
import sys

# Method to read the given species pair file and put them into a list of tuples
def pairs_of_eukaryotic_species(user_input_data='oma-pair_eukaryote.txt'):
    speciesPairs = []
    with open(user_input_data, 'r') as user_input_species_list:
        for line in user_input_species_list:
            new_line = line.split()
            speciesPairs.append((new_line[0], new_line[1]))
    return speciesPairs

# Method to read an amino acid alignment file and return a dictionary of sequence ID and sequence
def read_alignment(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()

    alignment_seq = {}
    seq_id = ""
    sequence = ""
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if seq_id != "":
                alignment_seq[seq_id] = sequence
            seq_id = line[1:].strip()
            sequence = ""
        else:
            sequence += line

    if seq_id != "":
        alignment_seq[seq_id] = sequence

    return alignment_seq

# Method to align DNA sequences to amino acid sequences
def align_dna_sequences(protein_alignment, dna_sequences):
    aligned_dna_sequences = {}

    for seq_id in protein_alignment:
        protein_seq = protein_alignment[seq_id]
        if seq_id in dna_sequences:
            dna_seq = dna_sequences[seq_id]

            aligned_dna_seq = ""
            index = 0
            for aa in protein_seq:
                if aa == "-":
                    aligned_dna_seq += "---"  # add "---" for each "-" in protein sequence
                else:
                    codon = dna_seq[index:index + 3]
                    aligned_dna_seq += codon
                    index += 3

            aligned_dna_sequences[seq_id] = aligned_dna_seq

    return aligned_dna_sequences

# Method to write a DNA sequence to the alignment file
def write_alignment(dna_id, dna_sequence, file_path):
    with open(file_path, "a") as outfile:
        outfile.write("> {}\n".format(dna_id))
        outfile.write("{}\n".format(dna_sequence))

# Method to align DNA sequences to amino acid sequences for a given species pair
def align_dna_to_aa(pair):
    sp1 = pair[0]
    sp2 = pair[1]

    print("Working on {} vs {}".format(sp1, sp2))

    # specify the folder containing the aligned protein files
    protein_folder = 'pairs_fasta/' + sp1 + '_' + sp2 + "/*_aa_*_aligned.fasta"
    # get the list of protein files
    protein_files = glob.glob(protein_folder)

    # specify the files with the raw DNA sequences
    species1_file = 'species_fasta/' + sp1 + '_dna' + '.fasta'
    species2_file = 'species_fasta/' + sp2 + '_dna' + '.fasta'

    # read the DNA sequences for each species from the DNA files
    dna_sequences = {}
    for dna_file in (species1_file, species2_file):
        dna_alignment = read_alignment(dna_file)
        dna_sequences.update(dna_alignment)

    for protein_file in tqdm.tqdm(protein_files):
        # read the aligned protein file
        protein_alignment = read_alignment(protein_file)

        # align the DNA sequences for the species pair
        aligned_dna_sequences = align_dna_sequences(protein_alignment, dna_sequences)

        # write the aligned DNA sequences to an output file
        output_file = protein_file.replace("_aa_", "_dna_")

        # clear the output file (in case it exists from a previous run)
        if os.path.exists(output_file):
            os.remove(output_file)

        for seq_id, dna_sequence in aligned_dna_sequences.items():
            write_alignment(seq_id, dna_sequence, output_file)

# Main function to orchestrate the whole process
def main():
    if len(sys.argv) == 1:
        # Step 1: Get the name of the file containing the species pairs
        userdata = input("Enter the file of species pairs: ")
    else:
        userdata = sys.argv[1]

    # Step 2: Read in the file
    species_pairs_list = pairs_of_eukaryotic_species(userdata)

    # Step 3: Align DNA to amino acid for each pair from user data
    for pair in species_pairs_list:
        align_dna_to_aa(pair)

if __name__ == "__main__":
    main()
