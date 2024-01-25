"""
primary calculation program

this program will read in a series of files containing a pair of aligned DNA sequences and calculate a LOT of things
from them, outputting the calculated values to a file, one output file per species pair containing one line per gene
pair
"""

import glob
import math
import sys


# this method reads in the given species pair file and puts them into a list of tuples
def pairs_of_eukaryotic_species(user_input_data='oma-pair_eukaryote.txt'):
    species_pairs = []
    with open(user_input_data, 'r') as user_input_species_list:
        for line in user_input_species_list:
            new_line = line.split()
            species_pairs.append((new_line[0], new_line[1]))
    return species_pairs


# def kappa_calc(total_nucleotides, a_count, t_count, g_count_average, c_count_average, purine_transitions,
#                pyrimidine_transitions, transversion_count):
def kappa_calc(total_nucleotides, A_frequency, T_frequency, G_frequency, C_frequency, purine_transitions,
               pyrimidine_transitions, transversion_count):
    # calculate each pi value/ frequency
    # A_frequency = a_count / total_nucleotides
    # T_frequency = t_count / total_nucleotides
    # G_frequency = g_count_average / total_nucleotides
    # C_frequency = c_count_average / total_nucleotides
    d, S, V = 0, 0, 0
    try:
        # purine and pyrimidine calculation
        R_frequency = A_frequency + G_frequency
        Y_frequency = C_frequency + T_frequency

        # k1 calculation below
        Kappa_1 = (2 * A_frequency * G_frequency) / R_frequency

        # k2 calculation below
        Kappa_2 = (2 * C_frequency * T_frequency) / Y_frequency

        # k3 calculation below
        Kappa_3 = 2 * ((R_frequency * Y_frequency) - ((A_frequency * G_frequency * Y_frequency) / R_frequency) - (
                (C_frequency * T_frequency * R_frequency) / Y_frequency))

        # proportion calculations
        R_proportion = purine_transitions / total_nucleotides
        Y_proportion = pyrimidine_transitions / total_nucleotides
        Q_proportion = transversion_count / total_nucleotides

        # w1 calculation below
        w_1 = 1 - (R_proportion / Kappa_1) - (Q_proportion / (2 * R_frequency))

        # w2 calculation below
        w_2 = 1 - (Y_proportion / Kappa_2) - (Q_proportion / (2 * Y_frequency))

        # w3 calculation below
        w_3 = 1 - (Q_proportion / (2 * R_frequency * Y_frequency))

        # S calculation below
        S = -1 * Kappa_1 * math.log(w_1) - Kappa_2 * math.log(w_2) - (Kappa_3 - 2 * R_frequency * Y_frequency) * math.log(w_3)
        S = max(0, S)  # value should not be negative
        # V calculation below
        V = -2 * R_frequency * Y_frequency * math.log(w_3)
        V = max(0, V)  # value should not be negative

        # distance
        d = -1 * Kappa_1 * math.log(w_1) - Kappa_2 * math.log(w_2) - Kappa_3 * math.log(w_3)
        d = max(0, d)  # value should not be negative

        # kappa approximation calculation below
        Kappa = (S / V) * ((R_frequency * Y_frequency) / (A_frequency * G_frequency + C_frequency * T_frequency))
        if S == 0:  # do not define Kappa if there are no transitions
            Kappa = -1
        return d, S, V, Kappa
    except (ZeroDivisionError, ValueError):
        if purine_transitions == 0 and pyrimidine_transitions == 0 and transversion_count == 0:  # sequences are identical
            return 0, 0, 0, -1
        elif transversion_count == 0:  # transitions but no transversions
            return d, S, V, -1
        return -1, -1, -1, -1


def filter_x(seq1, seq2):
    output1 = ""
    output2 = ""
    for i in range(len(seq1)):
        if (seq1[i] != "X") and (seq2[i] != "X"):
            output1 += seq1[i]
            output2 += seq2[i]
    return output1, output2


def filter_gaps(seq1, seq2):
    output1 = ""
    output2 = ""
    for i in range(len(seq1)):
        if (seq1[i] not in "-X") and (seq2[i] not in "-X"):
            output1 += seq1[i]
            output2 += seq2[i]
    return output1, output2


def filter_1st_position(seq1, seq2):
    output1 = seq1[::3]
    output2 = seq2[::3]
    output1, output2 = filter_gaps(output1, output2)
    return output1, output2


def filter_2nd_position(seq1, seq2):
    output1 = seq1[1::3]
    output2 = seq2[1::3]
    output1, output2 = filter_gaps(output1, output2)
    return output1, output2


def filter_3rd_position(seq1, seq2):
    output1 = seq1[2::3]
    output2 = seq2[2::3]
    output1, output2 = filter_gaps(output1, output2)
    return output1, output2


# define sets of synonymous codons for each amino acid
def is_4_fold_degenerate(codon):
    four_fold_degenerate_sets = [
        {"CTT", "CTC", "CTA", "CTG"},  # Leucine (L)
        {"GTT", "GTC", "GTA", "GTG"},  # Valine (V)
        {"TCT", "TCC", "TCA", "TCG"},  # Serine (S)
        {"CCT", "CCC", "CCA", "CCG"},  # Proline (P)
        {"ACT", "ACC", "ACA", "ACG"},  # Threonine (T)
        {"GCT", "GCC", "GCA", "GCG"},  # Alanine (A)
        {"CGT", "CGC", "CGA", "CGG"},  # Arginine (R)
        {"GGT", "GGC", "GGA", "GGG"},  # Glycine (G)
    ]

    # if a given codon is part of any of the 4-fold degenerate sets defined in four_fold_degenerate_sets,
    # returns True if it is, or False if it's not.
    for aa_set in four_fold_degenerate_sets:
        if codon in aa_set:
            return True
    return False


def filter_4_fold(seq1, seq2, filter_cpg: bool = False):
    output1 = ""  # initialize the output for seq1
    output2 = ""  # initialize the output for seq2

    for i in range(0, len(seq1), 3):
        codon1 = seq1[i:i+3]  # extract a codon of length 3 from seq1
        codon2 = seq2[i:i+3]  # extract a codon of length 3 from seq2
        # cpg checking
        check1a = seq1[i+1:i+3]  # check 2nd and 3rd position for a "CA" which could be deaminated CpG on reverse strand
        check2a = seq2[i+1:i+3]
        check1b = seq1[i+2:i+4]  # check 3rd pos and 1st of next codon for a "TG" which could be deaminated CpG
        check2b = seq2[i+2:i+4]
        # print(check1a, check1b, check2a, check2b)

        # check if either triplet is --- or if there is an X in it and if so skip
        if (codon1 != "---") and (codon2 != "---") and ("X" not in codon1) and ("X" not in codon2):
            # ask if codon1 and codon2 are part of the same 4-fold degenerate set
            if is_4_fold_degenerate(codon1) and is_4_fold_degenerate(codon2):
                if not filter_cpg or ((check1a != "CA") and (check1b != "TG") and
                                      (check2a != "CA") and (check2b != "TG")):
                    # if they are, output the 3rd position only
                    output1 += codon1[2]  # append the 3rd position of codon1 to the output for seq1
                    output2 += codon2[2]  # append the 3rd position of codon2 to the output for seq2

    return output1, output2  # return the filtered sequences containing the 3rd positions of 4-fold degenerate codons


def count_char_in_seqs(x, seq1, seq2):
    return seq1.count(x) + seq2.count(x)


def calculation_set_b(seq1, seq2):
    num_of_align_sites = len(seq1)

    # initialize counters
    identical_count = 0
    purine_transition_count = 0
    pyrimidine_transition_count = 0
    transversion_count = 0

    # iterate through each position in the alignment
    for i in range(num_of_align_sites):
        if seq1[i] == seq2[i]:
            # sequences are identical, increment identical count
            identical_count += 1
        elif (seq1[i] in "AG") and (seq2[i] in "AG"):  # purine transition
            purine_transition_count += 1
        elif (seq1[i] in "CT") and (seq2[i] in "CT"):  # pyrimidine transition
            pyrimidine_transition_count += 1
        else:
            # it must be a transversion, increment count
            transversion_count += 1

    return num_of_align_sites, identical_count, purine_transition_count, pyrimidine_transition_count, transversion_count


def add_b_results(prefix, data, n, id, purine, pyrimidine, transversion):
    data[prefix + " count sites"] = n
    data[prefix + " count identical"] = id
    data[prefix + " count purine transitions"] = purine
    data[prefix + " count pyrimidine transitions"] = pyrimidine
    data[prefix + " count transversions"] = transversion


def add_kappa_results(prefix, d, s, v, kappa, data):
    data[prefix + " TN distance"] = d
    data[prefix + " transition distance"] = s
    data[prefix + " transversion distance"] = v
    data[prefix + " kappa"] = kappa


def pair_calculations(pair):
    sp1 = pair[0]
    sp2 = pair[1]
    # get list of all DNA files for this pair
    filelist = glob.glob("pairs_fasta/" + sp1 + "_" + sp2 + "/*_dna_*_aligned.fasta")

    # open file for output
    with open("pairs_calcs/" + sp1 + "_" + sp2 + "_calcs.txt", "w") as outfile:
        # write a header line with column labels: we're not even sure what all of those will be yet
        headers = ["Species1",
                   "Species2",
                   "Sp1Gene",
                   "Sp2Gene",
                   "Ortholog Pair #",
                   "Aligned Length (w/gaps)",
                   "FreqA (total)",
                   "FreqC (total)",
                   "FreqG (total)",
                   "FreqT (total)"]

        repeat_headers = [" count sites",
                          " count identical",
                          " count purine transitions",
                          " count pyrimidine transitions",
                          " count transversions",
                          " TN distance",
                          " transition distance",
                          " transversion distance",
                          " kappa"]
        prefixes = ["ungapped", "4-fold", "4-fold cpg", "1st pos", "2nd pos", "3rd pos"]
        for p in prefixes:
            for h in repeat_headers:
                headers.append(p + h)
        # for h in repeat_headers:
        #     headers.append("ungapped" + h)
        # for h in repeat_headers:
        #     headers.append("1st pos" + h)
        # for h in repeat_headers:
        #     headers.append("2nd pos" + h)
        # for h in repeat_headers:
        #     headers.append("3rd pos" + h)

        outfile.write("\t".join(headers) + "\n")

        for filen in filelist:
            # use header labels as keys
            data = {"Species1": sp1, "Species2": sp2}

            # extract ortholog pair # of filen and add to dictionary
            filename_parts = filen.split("_")  # split the filename by underscores
            pair_number = filename_parts[5]  # extract the ortholog pair number from the last part of the filename before the extension
            data["Ortholog Pair #"] = pair_number

            # STEP: read filen and get the two sequences from it. They are already aligned DNA sequences
            with open(filen, "r") as dna_file:
                lines = dna_file.readlines()
                seq1_header = lines[0].strip()  # the first line contains the header for seq1
                seq1 = lines[1].strip()  # the second line contains the sequence for seq1
                seq2_header = lines[2].strip()  # the third line contains the header for seq2
                seq2 = lines[3].strip()  # the fourth line contains the sequence for seq2
                
            # extract the sequence name from the header
            seq1_name = seq1_header[1:] # header starts with ">"
            seq2_name = seq2_header[1:] # header starts with ">"
                
            data["Sp1Gene"] = seq1_name 
            data["Sp2Gene"] = seq2_name

            # STEP: calculations on entirety of sequence (except for sites with an X)

            filtered1, filtered2 = filter_x(seq1, seq2)
            data["Aligned Length (w/gaps)"] = len(filtered1)
            cnt_gaps = count_char_in_seqs("-", filtered1, filtered2)
            non_gap_chars = 2*len(filtered1) - cnt_gaps
            fa = count_char_in_seqs("A", filtered1, filtered2) / non_gap_chars
            fc = count_char_in_seqs("C", filtered1, filtered2) / non_gap_chars
            fg = count_char_in_seqs("G", filtered1, filtered2) / non_gap_chars
            ft = count_char_in_seqs("T", filtered1, filtered2) / non_gap_chars
            data["FreqA (total)"] = fa
            data["FreqC (total)"] = fc
            data["FreqG (total)"] = fg
            data["FreqT (total)"] = ft

            # create subseqs where any sites with gaps have been removed. From these subseqs calculate:
            filtered1, filtered2 = filter_gaps(seq1, seq2)
            n, id_cnt, purine_cnt, pyrimidine_cnt, transversion_cnt = calculation_set_b(filtered1, filtered2)
            add_b_results("ungapped", data, n, id_cnt, purine_cnt, pyrimidine_cnt, transversion_cnt)
            d, s, v, kappa = kappa_calc(n, fa, ft, fg, fc, purine_cnt, pyrimidine_cnt, transversion_cnt)
            add_kappa_results("ungapped", d, s, v, kappa, data)

            # create subseqs from original seqs that represent only 4-fold degenerate sites
            # from these calculate, again
            filtered1, filtered2 = filter_4_fold(seq1, seq2)
            n, id_cnt, purine_cnt, pyrimidine_cnt, transversion_cnt = calculation_set_b(filtered1, filtered2)
            add_b_results("4-fold", data, n, id_cnt, purine_cnt, pyrimidine_cnt, transversion_cnt)
            d, s, v, kappa = kappa_calc(n, fa, ft, fg, fc, purine_cnt, pyrimidine_cnt, transversion_cnt)
            add_kappa_results("4-fold", d, s, v, kappa, data)

            # create subseqs from original seqs that represent only 4-fold degenerate sites with potential CpG
            # sites removed, and then calculate, again

            filtered1, filtered2 = filter_4_fold(seq1, seq2, True)
            n, id_cnt, purine_cnt, pyrimidine_cnt, transversion_cnt = calculation_set_b(filtered1, filtered2)
            add_b_results("4-fold cpg", data, n, id_cnt, purine_cnt, pyrimidine_cnt, transversion_cnt)
            d, s, v, kappa = kappa_calc(n, fa, ft, fg, fc, purine_cnt, pyrimidine_cnt, transversion_cnt)
            add_kappa_results("4-fold cpg", d, s, v, kappa, data)

            # MAYBE
            # create subseqs from original seqs that represent only 2-fold degenerate sites
            # from these calculate, again

            #   CALC SET B
            #   CALC SET C
            #   CALC SET D

            # MAYBE
            # create subseqs from original seqs that represent only 0-fold degenerate sites
            # from these calculate, again

            #   CALC SET B
            #   CALC SET C
            #   CALC SET D



            # MAYBE
            # create subseqs from original seqs that represent only 1st codon positions
            filtered1, filtered2 = filter_1st_position(seq1, seq2)
            n, id_cnt, purine_cnt, pyrimidine_cnt, transversion_cnt = calculation_set_b(filtered1, filtered2)
            add_b_results("1st pos", data, n, id_cnt, purine_cnt, pyrimidine_cnt, transversion_cnt)
            d, s, v, kappa = kappa_calc(n, fa, ft, fg, fc, purine_cnt, pyrimidine_cnt, transversion_cnt)
            add_kappa_results("1st pos", d, s, v, kappa, data)

            # MAYBE
            # create subseqs from original seqs that represent only 2nd codon positions
            filtered1, filtered2 = filter_2nd_position(seq1, seq2)
            n, id_cnt, purine_cnt, pyrimidine_cnt, transversion_cnt = calculation_set_b(filtered1, filtered2)
            add_b_results("2nd pos", data, n, id_cnt, purine_cnt, pyrimidine_cnt, transversion_cnt)
            d, s, v, kappa = kappa_calc(n, fa, ft, fg, fc, purine_cnt, pyrimidine_cnt, transversion_cnt)
            add_kappa_results("2nd pos", d, s, v, kappa, data)

            # MAYBE
            # create subseqs from original seqs that represent only 3rd codon positions
            filtered1, filtered2 = filter_3rd_position(seq1, seq2)
            n, id_cnt, purine_cnt, pyrimidine_cnt, transversion_cnt = calculation_set_b(filtered1, filtered2)
            add_b_results("3rd pos", data, n, id_cnt, purine_cnt, pyrimidine_cnt, transversion_cnt)
            d, s, v, kappa = kappa_calc(n, fa, ft, fg, fc, purine_cnt, pyrimidine_cnt, transversion_cnt)
            add_kappa_results("3rd pos", d, s, v, kappa, data)

            #  write a single line to the output file with ALL calculated values for this sequence pair,
            #     essentially everything in the dictionary, matching the order of the headers above

            # create a list with strings containing all output values for the line
            #   the first five entries in the table will be strings (I think)
            outline = [data[x] for x in headers[:5]]
            #   everything else will be a number. Some will be integers and some real; the function num2str attempts
            #   to optimally convert based on which is which
            outline.extend([num2str(data[x]) for x in headers[5:]])
            outfile.write("\t".join(outline) + "\n")


def num2str(x) -> str:
    """
    A little function to convert a number to a string, treating integers and floating-point numbers differently
    """
    if isinstance(x, int):
        return str(x)
    else:
        return format(x, "0.6f")  # 6 decimal places should be adequate


def main():
    if len(sys.argv) == 1:
        # step 1 (Get the name of the file which contains the species pairs you want to analyze)
        userdata = input("Enter the file of species pairs: ")
    else:
        userdata = sys.argv[1]

    # step 2 (reading in the file)
    species_pairs_list = pairs_of_eukaryotic_species(userdata)

    # step 3 (perform calculations for each pair of species)
    for pair in species_pairs_list:
        pair_calculations(pair)


if __name__ == "__main__":
    main()
