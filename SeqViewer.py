import re
import tkinter as tk
from tkinter import filedialog, Text, Listbox, Scrollbar, Toplevel, Label, Button
from collections import OrderedDict
import os
from tkinter.ttk import Treeview
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import regex
from NeedlemanWunsch_class_updated import Needleman_Wunsch
import cgi
import mysql
import mysql.connector
from mysql.connector import Error
import tkinter.messagebox as messagebox



host = "localhost"
user = "root"
password = "&%Bn96=mdQe4"
database = "liuac"



def is_empty_file(filepath):
    return os.path.getsize(filepath) == 0
def sort_treeview(tree, col, descending):
    if col == "GC Content":
        data = [(float(tree.set(item, col).rstrip("%")), item) for item in tree.get_children('')]
    else:
        data = [(tree.set(item, col).lower() if col not in ["Length", "A", "T", "G", "C"] else int(tree.set(item, col)), item) for item in tree.get_children('')]

    data.sort(reverse=descending)
    for index, (val, item) in enumerate(data):
        tree.move(item, '', index)
    tree.heading(col, command=lambda: sort_treeview(tree, col, not descending))



def print_with_ruler(seq_name, description, sequence, NucleotidesPerLine, spacer):
    """Prints sequence information with a ruler.

    Args:
        seq_name (str): The name of the sequence.
        description (str): The description or additional information about the sequence.
        sequence (str): The sequence bases.
        NucleotidesPerLine (int): The maximum number of nucleotides to print per line.
        spacer (bool): A flag indicating whether to add spaces between nucleotides.
    """
    result = ""
    if (
        spacer == False
    ):  # If spacer is set to False, do not add spaces between nucleotides.
        print(
            # ">" + seq_name, description, "\n"
        )  # Print the sequence name and description.
        result = f">{seq_name} {description}\n"
        repeat_count = NucleotidesPerLine // 10

        # Generate a string of repeated digits from 1 to 10 (0) based on repeat_count.
        repeated_string = "".join(["1234567890"] * repeat_count)

        # Create a line header with numbers and spaces.
        line_header = " ".join(["Line", repeated_string])

        # Print the header row with line numbers and ruler.
        # print(f"{1 :> 15}", end="")  # Print the first line number.
        result += f"{1 :> 15}"
        for k in range(repeat_count - 1):
            # print(f"{k + 2 :> 10}", end="")  # Print subsequent line numbers.
            result += f"{k + 2 :> 10}"
        # print()  # Move to the next line.
        # print(line_header)  # Print the ruler.
        result += f"\n{line_header}"

        # Print the sequence with line numbers and without spaces.
        for i in range(0, len(sequence), NucleotidesPerLine):
            chunk = sequence[i : i + NucleotidesPerLine]

            # Remove spaces between every 10 characters in the chunk.
            chunk_without_spaces = "".join(
                [chunk[j : j + 10] for j in range(0, len(chunk), 10)]
            )

            # Print the line number and sequence chunk.
            # print(f"{i // NucleotidesPerLine + 1 :> 4}", chunk_without_spaces)
            result += f"\n{i // NucleotidesPerLine + 1 :> 4} {chunk_without_spaces}"
    else:
        # print(
        #     ">" + seq_name, description, "\n"
        # )  # Print the sequence name and description.
        result = f">{seq_name} {description}\n"
        repeat_count = NucleotidesPerLine // 10

        # Generate a string of repeated digits from 1 to 10 (0 based on repeat_count.
        repeated_string = " ".join(["1234567890"] * repeat_count)

        # Create a line header with numbers and spaces.
        line_header = " ".join(["Line", repeated_string])

        # Print the header row with line numbers and ruler.
        # print(f"{1 :> 15}", end="")  # Print the first line number.
        result += f"{1 :> 15}"
        for k in range(repeat_count - 1):
            # print(
            #     f"{k + 2 :> 11}", end=""
            # )  # Print subsequent line numbers with extra spacing.
            result += f"{k + 2 :> 11}"
        # print()  # Move to the next line.
        # print(line_header)  # Print the ruler.
        result += f"\n{line_header}"

        # Print the sequence with line numbers and spaces.
        for i in range(0, len(sequence), NucleotidesPerLine):
            chunk = sequence[i : i + NucleotidesPerLine]

            # Add spaces between every 10 characters in the chunk.
            chunk_with_spaces = " ".join(
                [chunk[j : j + 10] for j in range(0, len(chunk), 10)]
            )

            # Print the line number and sequence chunk.
            # print(f"{i // NucleotidesPerLine + 1 :> 4}", chunk_with_spaces)
            result += f"\n{i // NucleotidesPerLine + 1: >4} {chunk_with_spaces}"
    return result
def detect_homopolymer(sequence):
    """Detect homopolymer sequences in a DNA sequence.

    Args:
        sequence (str): The input DNA sequence to search for homopolymers.

    Returns:
        list: A list of homopolymer sequences in the format 'start-end_length_nucleotide'.
    """
    homopolymers = []  # List to store detected homopolymers

    for nuc in ['A', 'T', 'G', 'C']:        
        pattern = r"({nuc}{3,}{nuc}{s<=1}{nuc}{3,})"
        pattern = pattern.replace("{nuc}", nuc)  # Replace the placeholder with the nucleotide

        regex_pattern = regex.compile(pattern)
        for match in regex_pattern.finditer(sequence):
            start = match.start()
            end = match.end() - 1  # Adjust the end position
            length = end - start + 1
            match_count = match.group(0).count(nuc)  # Count the occurrences of the nucleotide in the match
            if match_count >= 10:  # Check if the homopolymer length is at least 10
                match_count = match.group(0).count(nuc)  # Count the occurrences of the nucleotide in the match
                # print(match_count)
                homopolymers.append(f'{start}-{end}_{length}_{nuc}')

    return homopolymers

def motif_search(sequence, target):
    """Search for a target motif within a sequence and return the positions of detected motifs.

    Args:
        sequence (str): The input sequence to search for motifs.
        target (str): The target motif to search for within the sequence.

    Returns:
        list: A list of motif positions in the format 'start-end_length' for each detected motif.
    """
    detected_motifs = []  # List to store detected motif positions
    seq_len = len(sequence)  # Length of the input sequence
    target_len = len(target)  # Length of the target motif
    i = 0  # Initialize the index for iterating through the sequence

    while i < seq_len:
        if sequence[i:i + target_len] == target:
            # If the current substring matches the target motif, record its position
            start_pos = i
            end_pos = i + target_len - 1
            motif_length = target_len
            detected_motifs.append(f"{start_pos}-{end_pos}_{motif_length}")
            i += target_len  # Skip the length of the motif to avoid overlapping detections
        else:
            i += 1  # Move to the next position in the sequence if no match is found

    return detected_motifs  # Return the list of detected motif positions

def modify_sequence(sequence, homopolymers, motifs):
    """Modify a given DNA sequence by converting specified regions to lowercase letters.

    Args:
        sequence (str): The input DNA sequence to be modified.
        homopolymers (list): A list of homopolymer information (start-end_length_letter) to be converted to lowercase.
        motifs (list): A list of motif information (start-end_length) to be converted to lowercase.

    Returns:
        str: The modified DNA sequence with specified regions in lowercase.
    """
    modified_sequence = list(sequence)  # Convert the sequence to a list of characters for easy modification

    if motifs is not None:
        for motif_info in motifs:
            start_end, length = motif_info.split('_')
            start, end = map(int, start_end.split('-'))

            for i in range(start, end + 1):
                modified_sequence[i] = sequence[i].lower()

    if homopolymers is not None:
        for homopolymer_info in homopolymers:
            start_end, length, letter = homopolymer_info.split('_')
            start, end = map(int, start_end.split('-'))

            for i in range(start, end + 1):
                modified_sequence[i] = sequence[i].lower()

    return ''.join(modified_sequence)  # Convert the modified list back to a string

def print_with_ruler2(sequence, NucleotidesPerLine, spacer):
    """Prints sequence information with a ruler.

    Args:
        seq_name (str): The name of the sequence.
        description (str): The description or additional information about the sequence.
        sequence (str): The sequence bases.
        NucleotidesPerLine (int): The maximum number of nucleotides to print per line.
        spacer (bool): A flag indicating whether to add spaces between nucleotides.
    """
    if (
        spacer == False
    ):  # If spacer is set to False, do not add spaces between nucleotides.
        repeat_count = NucleotidesPerLine // 10

        # Generate a string of repeated digits from 1 to 10 (0) based on repeat_count.
        repeated_string = "".join(["1234567890"] * repeat_count)

        # Create a line header with numbers and spaces.
        line_header = " ".join(["Line", repeated_string])

        # Print the header row with line numbers and ruler.
        print(f"{1 :> 15}", end="")  # Print the first line number.
        for k in range(repeat_count - 1):
            print(f"{k + 2 :> 10}", end="")  # Print subsequent line numbers.
        print()  # Move to the next line.
        print(line_header)  # Print the ruler.

        # Print the sequence with line numbers and without spaces.
        for i in range(0, len(sequence), NucleotidesPerLine):
            chunk = sequence[i : i + NucleotidesPerLine]

            # Remove spaces between every 10 characters in the chunk.
            chunk_without_spaces = "".join(
                [chunk[j : j + 10] for j in range(0, len(chunk), 10)]
            )

            # Print the line number and sequence chunk.
            print(f"{i // NucleotidesPerLine + 1 :> 4}", chunk_without_spaces)
    else:
        repeat_count = NucleotidesPerLine // 10

        # Generate a string of repeated digits from 1 to 10 (0 based on repeat_count.
        repeated_string = " ".join(["1234567890"] * repeat_count)

        # Create a line header with numbers and spaces.
        line_header = " ".join(["Line", repeated_string])

        # Print the header row with line numbers and ruler.
        print(f"{1 :> 15}", end="")  # Print the first line number.
        for k in range(repeat_count - 1):
            print(
                f"{k + 2 :> 11}", end=""
            )  # Print subsequent line numbers with extra spacing.
        print()  # Move to the next line.
        print(line_header)  # Print the ruler.

        # Print the sequence with line numbers and spaces.
        for i in range(0, len(sequence), NucleotidesPerLine):
            chunk = sequence[i : i + NucleotidesPerLine]

            # Add spaces between every 10 characters in the chunk.
            chunk_with_spaces = " ".join(
                [chunk[j : j + 10] for j in range(0, len(chunk), 10)]
            )

            # Print the line number and sequence chunk.
            print(f"{i // NucleotidesPerLine + 1 :> 4}", chunk_with_spaces)


def print_targets(sequence, list_homo, list_motif):
    """Print detected homopolymers and motifs in a modified DNA sequence, along with the modified sequence.

    Args:
        sequence (str): The original DNA sequence.
        list_homo (list): A list of detected homopolymers.
        list_motif (list): A list of detected motifs.
    """
    new_sequence = modify_sequence(sequence, list_homo, list_motif)

    if list_homo is not None:
        print(f"Homopolymer Search:", end=" ")
        for idx, homopolymer in enumerate(list_homo, start=1):
            print(f"{idx}={homopolymer}", end=" ")
        print()

    if list_motif is not None:
        print(f"Motif Search for {config_dict['motifSearchTarget']}:", end=" ")
        for idx, motif in enumerate(list_motif, start=1):
            print(f"{idx}={motif}", end=" ")
        print()

    print()
    print_with_ruler2(new_sequence, int(config_dict["NucleotidesPerLine[50|100]"]), config_dict["DoYouNeedSpaceSeperator[Y|N]"] == "Y")
    print()

def cpg_island(sequence):
    """Identify CpG islands in a given DNA sequence.

    Args:
        sequence (str): The input DNA sequence.

    Returns:
        dict: A dictionary containing information about identified CpG islands.
            The keys are unique identifiers for each island (1-based index), and
            the values are strings in the format "start-end_length" where:
                - start: The start position of the CpG island.
                - end: The end position of the CpG island (exclusive).
                - length: The length of the CpG island.
    """
    i = 0
    count = 0
    cpg_dict = {}

    # Iterate through the sequence, considering pairs of consecutive nucleotides.
    while i in range(len(sequence) - 1):
        dinucleotide = sequence[i : i + 2]

        # Check if the current dinucleotide is "CG" (C followed by G).
        if dinucleotide == "CG":
            start = i  # Store the start position of the CpG island.

            # Continue checking for consecutive "CG" dinucleotides.
            while dinucleotide == "CG":
                dinucleotide = sequence[i : i + 2]
                i = i + 2  # Move to the next pair of nucleotides within the island.

            end = i - 2  # Store the end position of the CpG island.

            # Check if the island has a length of 6 or more.
            if end - start >= 6:
                count = count + 1

                # Create a range string (e.g., "start-end") and a length string.
                ranges = "-".join([str(start), str(end - 1)])
                value = "_".join([ranges, str(end - start)])

                # Store the CpG island information in the cpg_dict dictionary.
                cpg_dict[count] = value

        i = i + 2  # Move to the next pair of nucleotides in the sequence.

    # Return a dictionary containing information about identified CpG islands.
    return cpg_dict

def translation(dna_sequence):
    """Translate a DNA sequence into a protein sequence.

    Args:
        dna_sequence (str): The input DNA sequence to be translated.

    Returns:
        str: The resulting protein sequence.
    """
    trans_dic = {
        "UUU": "F",
        "UUC": "F",
        "UUA": "L",
        "UUG": "L",
        "CUU": "L",
        "CUC": "L",
        "CUA": "L",
        "CUG": "L",
        "AUU": "I",
        "AUC": "I",
        "AUA": "I",
        "AUG": "M",
        "GUU": "V",
        "GUC": "V",
        "GUA": "V",
        "GUG": "V",
        "UCU": "S",
        "UCC": "S",
        "UCA": "S",
        "UCG": "S",
        "CCU": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "ACU": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "GCU": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "UAU": "Y",
        "UAC": "Y",
        "UAA": "*",
        "UAG": "*",
        "CAU": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "AAU": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "GAU": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "UGU": "C",
        "UGC": "C",
        "UGA": "*",
        "UGG": "W",
        "CGU": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "AGU": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GGU": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }
    rna_sequence = dna_sequence.replace("T", "U")
    protein = []
    result = ""
    for i in range(0, len(rna_sequence) - 2, 3):
        protein.append(trans_dic[rna_sequence[i : i + 3]])
    result = result.join(protein)
    return result

def reverse_complement(sequence):
    """Calculate the reverse complementary DNA strand.

    Args:
        sequence (str): The input DNA sequence.

    Returns:
        str: The reverse complementary DNA strand.
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    reverse_seq = sequence[::-1]
    reverse_comp_seq = "".join(complement[base] for base in reverse_seq)
    return reverse_comp_seq

def translation_6_frame(dna_sequence):
    """
    Translate a DNA sequence in all six reading frames.

    Args:
        dna_sequence (str): The input DNA sequence to be translated.

    Returns:
        Tuple of six protein sequences (str): The protein sequences translated from the six reading frames.
    """
    # Translate the original DNA sequence in the first reading frame
    protein_1 = translation(dna_sequence)

    # Translate the DNA sequence shifted by one base to the right (second reading frame)
    protein_2 = translation(dna_sequence[1:])

    # Translate the DNA sequence shifted by two bases to the right (third reading frame)
    protein_3 = translation(dna_sequence[2:])

    # Calculate the reverse complement of the DNA sequence
    rev_comp = reverse_complement(dna_sequence)
    # Translate the reverse complement in the first reading frame
    protein_4 = translation(rev_comp)

    # Translate the reverse complement shifted by one base to the right (fifth reading frame)
    protein_5 = translation(rev_comp[1:])

    # Translate the reverse complement shifted by two bases to the right (sixth reading frame)
    protein_6 = translation(rev_comp[2:])

    # Return a tuple containing the six protein sequences
    return (protein_1, protein_2, protein_3, protein_6[::-1], protein_5[::-1], protein_4[::-1])


# Need to revise printSeqFragment, one of the fragments on the top
# or bottom is incorrect
def print_seq_fragment(seq_fragment, start, end):
    """
    Print a sequence fragment with annotations and translations.

    Args:
        seq_fragment (str): The sequence fragment to be printed.
        start (int): The start index of the fragment in the original sequence.
        end (int): The end index of the fragment in the original sequence.
    """
    
    # Calculate the length of the sequence fragment
    seq_length = end - start + 1

    # Extract the sequence from the fragment
    sequence = seq_fragment[start - 1 : end]

    # Translate the sequence in all six reading frames
    proteins = translation_6_frame(sequence)

    # Print the translations in the first three reading frames
    for p in proteins[3:]:
        print(format_seq_frag(sequence, p))

    # Print the reverse complement of the sequence
    print(reverse_complement(sequence)[::-1], end="\n")

    # Print a line of vertical bars to separate the annotations
    for i in range(seq_length):
        print("|", end="")
    print()

    # Print a line with sequence range annotation
    print(f"<{start}{'-' * (end - start - len(str(start)) - len(str(end)) - 1)}{end}>")

    # Print another line of vertical bars
    for i in range(seq_length):
        print("|", end="")
    print()

    # Print the original sequence
    print(sequence)
    
    # Print the translations in the last three reading frames
    for p in proteins[0:3]:
        print(format_seq_frag(sequence, p))
    print()

    
def alignment(seq_name_list, seq_list):
    """Perform sequence alignment and print the results.

    Args:
        seq_name_list (list): A list of sequence names.
        seq_list (list): A list of sequences to align.
    """
    # Ensure there are at least two sequences for alignment
    if len(seq_list) < 2:
        print("At least two sequences are required for alignment.")
        return

    # Define the NeedlemanWunsch object for sequence alignment
    # Align the first sequence with all other sequences
    first_sequence = seq_list[0]
    for i in range(1, len(seq_list)):
        aligner = Needleman_Wunsch(first_sequence, seq_list[i])
        aligned_seq1, aligned_seq2 = aligner.give_final_result()
        aligned_seq1 = seq_name_list[0] + "\t" + aligned_seq1
        aligned_seq2 = seq_name_list[i] + "\t" + aligned_seq2
        # Call the process_aligned_seq function to print alignments and CIGAR string
        process_aligned_seq(aligned_seq1, aligned_seq2)


def process_aligned_seq(aligned_seq1, aligned_seq2):
    """Process and print alignment information between two aligned sequences.

    Args:
        aligned_seq1 (str): First aligned sequence with a sequence name.
        aligned_seq2 (str): Second aligned sequence with a sequence name.
    """
    def is_pyrimidine(char):
        return char in ('C', 'T')

    def is_purine(char):
        return char in ('A', 'G')
    
    print(aligned_seq1)
    # Initialize variables for match, mismatch, insertion, and deletion counts
    match_count = 0
    mismatch_count = 0
    insertion_count = 0
    deletion_count = 0
    pyr_pur_count = 0

    alignment_line = ""  # Initialize the middle line of the alignment
    cigar_string = ""    # Initialize the CIGAR string

    # Calculate match, mismatch, insertion, and deletion counts and build the CIGAR string
    current_operation = ""  # Track the current operation (M, I, D)
    current_count = 0  # Track the current count of contiguous operations
    
    for char1, char2 in zip(aligned_seq1.split("\t")[1], aligned_seq2.split("\t")[1]):
        if char1 == char2:
            match_count += 1
            alignment_line += "|"  # Match symbol
            if current_operation != "M":
                if current_operation:
                    cigar_string += str(current_count) + current_operation
                current_operation = "M"
                current_count = 1
            else:
                current_count += 1
        elif char1 == '-':
            insertion_count += 1
            alignment_line += " "  # Insertion symbol
            if current_operation != "I":
                if current_operation:
                    cigar_string += str(current_count) + current_operation
                current_operation = "I"
                current_count = 1
            else:
                current_count += 1
        elif char2 == '-':
            deletion_count += 1
            alignment_line += " "  # Deletion symbol
            if current_operation != "D":
                if current_operation:
                    cigar_string += str(current_count) + current_operation
                current_operation = "D"
                current_count = 1
            else:
                current_count += 1
        else:
            if is_pyrimidine(char1) and is_pyrimidine(char2):
                alignment_line += ":"  # Colon for pyrimidines
                pyr_pur_count += 1
            elif is_purine(char1) and is_purine(char2):
                alignment_line += ":"  # Colon for purines
                pyr_pur_count += 1
            else:
                alignment_line += " "  # Mismatch symbol
            mismatch_count += 1
            if current_operation != "M":
                if current_operation:
                    cigar_string += str(current_count) + current_operation
                current_operation = "M"
                current_count = 1
            else:
                current_count += 1

    # Append the last operation to the CIGAR string
    if current_operation:
        cigar_string += str(current_count) + current_operation

    # Calculate match, mismatch, insertion, and deletion ratios
    total_length = len(aligned_seq1)
    match_ratio = match_count / total_length

    alignment_spaces = len(aligned_seq1.split("\t")[0]) * ' '
    alignment_line = alignment_spaces + "\t" + alignment_line
    print(alignment_line)
    print(aligned_seq2, end="\n\n")
    print(f"Identities={match_count}/{total_length} ({match_ratio * 100:.1f}%)", end= " ")
    print(f"Similarity={(match_count + pyr_pur_count)}/{total_length} ({(match_count + pyr_pur_count) / total_length * 100:.1f}%)", end=" ")
    print(f"Gaps={insertion_count + deletion_count}/{total_length} ({(insertion_count + deletion_count) / total_length * 100:.1f}%)", end=" ")
    print(f"CIGAR = {cigar_string}", end="\n\n\n")


def positional_matrix(seq_list):
    """Generate a positional matrix to represent the distribution of nucleotides at each position in the sequences.

    Args:
        seq_list (list): A list of DNA sequences.

    Returns:
        numpy.ndarray: A positional matrix where each row represents a nucleotide (A, T, G, C)
                      and each column represents a position in the sequences.
    """
    # Get the length of the sequences
    seq_length = len(seq_list[0])
    nucleotides = ['A', 'T', 'G', 'C']
    nucleotides_size = len(nucleotides)

    # Initialize an empty positional matrix filled with zeros
    pos_matrix = np.zeros((nucleotides_size, seq_length), dtype=int)

    # Create a mapping from characters to row indices in the matrix
    char_to_index = {'A': 0, 'T': 1, 'G': 2, 'C': 3}

    # Iterate through the sequences character by character
    for sequence in seq_list:
        for position, char in enumerate(sequence):
            if char in char_to_index:
                row_index = char_to_index[char]
                pos_matrix[row_index, position] += 1

    return pos_matrix

def show_positional_matrix(my_array):
    """Display and visualize a positional matrix, and save it as a stacked bar chart image.

    Args:
        my_array (numpy.ndarray): A positional matrix representing the distribution of nucleotides at each position.
    """
    # Get the number of positions from the positional matrix
    nucleotide_size, num_positions = my_array.shape
    left_align = len(str(num_positions))

    # Print position headers
    print("Position", end=" ")
    for i in range(1, num_positions + 1):
        print("{: <{width}}".format(i, width=left_align), end=" ")
    print()

    # Print a line separator
    print("--------", end=" ")
    for i in range(num_positions):
        print(f"{left_align * '-'}", end=" ")
    print()

    # Print the character counts for each nucleotide
    nucleotides = ['A', 'T', 'G', 'C']
    for row_index, nucleotide in enumerate(nucleotides):
        print(nucleotide, end=f"{8 * ' '}")
        for position in range(num_positions):
            print(f"{my_array[row_index, position]}", end=f"{left_align * ' '}")
        print()
        
        
    
    # Create a color map for nucleotides
    colors = {'A': 'blue', 'T': 'red', 'G': 'gray', 'C': 'orange'}
    
    # Calculate nucleotide percentages
    total_counts = np.sum(my_array, axis=0)
    percentages = (my_array / total_counts) * 100
    plt.figure(figsize=(16, 9))
    
    # Plot a stacked bar graph
    positions = range(1, num_positions + 1)
    nucleotides = ['A', 'T', 'G', 'C']
    bottom = np.zeros(num_positions)  # Initialize the bottom positions for stacking

    for nucleotide in nucleotides:
        counts = percentages[nucleotides.index(nucleotide)]
        plt.bar(positions, counts, label=nucleotide, color=colors[nucleotide], bottom=bottom)
        bottom += counts

    plt.title('Positional Nucleotide Matrix')
    
    plt.yticks(range(0, 101, 10))

    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=4, frameon=False, borderaxespad=0.5)
    
    # Save the stacked bar chart as a PNG file
    plt.savefig('PNM.jpg', format='JPG')
  
def codon_profile(sequence):
    """Generate a codon profile for a given DNA sequence.

    Args:
        sequence (str): A string representing a DNA sequence.

    Returns:
        dict: A dictionary where keys are codons (3-letter DNA sequences) and values are initialized to 0.
    """
    bases = ["A", "T", "C", "G"]
    codon_dict = {}  # Initialize an empty dictionary to store the codon profile.
    for base1 in bases:
        for base2 in bases:
            for base3 in bases:
                codon = (
                    base1 + base2 + base3
                )  # Generate a codon by concatenating three bases.
                codon_dict[codon] = 0  # Initialize the codon count to 0.

    for i in range(len(sequence)):
        for codon in codon_dict:
            if sequence[i : i + 3] == codon:
                codon_dict[codon] = codon_dict[codon] + 1
    return codon_dict

def codon_profile_print(codon_dict):
    """Print the codon profile based on the provided codon dictionary.

    Args:
        codon_dict (dict): A dictionary containing codon counts, where the keys
            are codons (e.g., "TAA", "GCC") and the values are their respective counts.
    """
    
    print("Codon Profile:", end="\n\n")
    print(f"{'2nd':>20}")
    print(f"{'-'*35:>43}")
    print("1st\tT\tC\tA\tG\t3rd")
    result = f"Codon Profile:\n\n"
    result += f"{'2nd':>20}\n{'-'*35:>43}\n1st\tT\tC\tA\tG\t3rd\n"
    # Define the order of bases
    bases = ["T", "C", "A", "G"]

    # Iterate over the first base
    for first_base in bases:
        row = first_base + "\t"

        # Counter to keep track of codons printed in the current line
        codons_printed = 0

        # Iterate over the second base
        for second_base in bases:
            # Iterate over the third base
            for third_base in bases:
                codon = first_base + third_base + second_base

                # Check if the codon exists in the dictionary
                if codon in codon_dict:
                    count = codon_dict[codon]

                # Append the codon and count to the row
                row += f"{codon}={count}\t"

                codons_printed += 1

                # Check if we have printed 4 codons, then start a new line
                if codons_printed == 4:
                    row += second_base  # Add the second base as 3rd column
                    print(row)
                    result += row + "\n"
                    row = "\t"  # Start a new row
                    codons_printed = 0
        print()
        result += "\n"
    return result
class ModifiedFASTAViewer:

    def upload_to_database(self):
        connection = create_connection(host, user, password, database)
        # Get data to insert into the database
        # name, description, sequence, length = self.get_header_sequence()

        # Establish a database connection
        try:
            with connection.cursor() as cursor:
                # SQL query to insert data into the SEQUENCE table
                sql = "INSERT INTO SEQUENCE (SEQUENCE_NAME, SEQUENCE_DESCRIPTION, SEQUENCE) VALUES (%s, %s, %s)"
                
                # Execute the query
                for kvp in self.sequences.items():
                
                    headers = kvp[0].split(">")
                
                    name = headers[1].split(" ")[0]
                
                    description = " ".join(headers[1].split(" ")[1:])
                    
                    sequence = kvp[1]
                    print(name, description, kvp[1])
                    cursor.execute(sql, (name, description, sequence))

            # Commit changes to the database
            connection.commit()

            # Show a success message box in the GUI
            success_message = "Data uploaded successfully!"
            messagebox.showinfo("Success", success_message)
            
        except Exception as e:
            print(f"Error: {e}")

        finally:
            # Close the database connection
            connection.close()

    def download_from_database(self):
        # Clear previous items in the Treeview
        for item in self.tree.get_children():
            self.tree.delete(item)

        # Establish a database connection
        connection = create_connection(host, user, password, database)
        
        try:
            with connection.cursor() as cursor:
                # SQL query to retrieve data from the SEQUENCE table
                sql = "SELECT SEQUENCE_NAME, SEQUENCE_DESCRIPTION, SEQUENCE FROM SEQUENCE"
                
                # Execute the query
                cursor.execute(sql)
                
                # Fetch all the results
                results = cursor.fetchall()

                # Populate the Treeview with sequence information
                for row in results:
                    name = row[0]
                    description = row[1]
                    sequence = row[2]
                    length = len(sequence)
                    a_count = sequence.count('A')
                    t_count = sequence.count('T')
                    g_count = sequence.count('G')
                    c_count = sequence.count('C')
                    gc_content = (g_count + c_count) / length * 100 if length > 0 else 0

                    self.tree.insert("", "end", values=(name, description, length, a_count, t_count, g_count, c_count, f"{gc_content:.2f}%"))
                    
                    # Update the sequences dictionary
                    header = f">{name} {description}"
                    self.sequences[header] = sequence

        except Exception as e:
            print(f"Error: {e}")

        finally:
            # Close the database connection
            connection.close()

        return
        
    def __init__(self, master):
        self.master = master
        self.master.title("FASTA Viewer")
        
        # Create a frame to hold the buttons
        self.button_frame = tk.Frame(self.master)
        self.button_frame.pack(pady=10)

        # Create a button to upload the FASTA file
        self.upload_button = Button(self.button_frame, text="Upload FASTA", command=self.upload_file)
        self.upload_button.pack(side="left", padx=5)

        # Create a button to upload to the database
        self.upload_to_database_button = Button(self.button_frame, text="Upload to Database", command=self.upload_to_database)
        self.upload_to_database_button.pack(side="left", padx=5)

        # Create a button to download from the database
        self.download_from_database_button = Button(self.button_frame, text="Download from Database", command=self.download_from_database)
        self.download_from_database_button.pack(side="left", padx=5)
        
        # List containing checkboxes
        self.cbvars = {}
        # Create a frame for checkboxes
        checkbox_frame = tk.Frame(self.master)
        checkbox_frame.pack()

        


        # Homopolymer checkbox
        hp_var = tk.IntVar()
        self.homopolymer = tk.Checkbutton(checkbox_frame, text="Homopolymer", onvalue=1, offvalue=0, variable=hp_var)
        self.homopolymer.grid(row=0, column=0, padx=5)
        self.cbvars["homopolymer"] = hp_var

        # CpG Island checkbox
        cpg_var = tk.IntVar()
        self.cpg_island = tk.Checkbutton(checkbox_frame, text="CpG Island", onvalue=1, offvalue=0, variable=cpg_var)
        self.cpg_island.grid(row=0, column=1, padx=5)
        self.cbvars["cpg"] = cpg_var

        # Motif Search checkbox
        motif_var = tk.IntVar()
        self.motif_search = tk.Checkbutton(checkbox_frame, text="Motif Search", onvalue=1, offvalue=0, variable=motif_var, command=self.toggle_motif_entry)
        self.motif_search.grid(row=0, column=2, padx=5)
        self.cbvars["motif"] = motif_var
        
        # Entry for motif search
        self.motif_entry = tk.Entry(checkbox_frame)
        self.motif_entry.grid(row=0, column=3, padx=5)
        self.motif_entry.grid_remove()  # Initially hide the entry

        # Codon Profile checkbox
        codon_var = tk.IntVar()
        self.codon_profile = tk.Checkbutton(checkbox_frame, text="Codon Profile", onvalue=1, offvalue=0, variable=codon_var)
        self.codon_profile.grid(row=0, column=4, padx=5)
        self.cbvars["codon"] = codon_var

        # Print Targets checkbox
        # ptargets_var = tk.IntVar()
        # self.printTargets = tk.Checkbutton(checkbox_frame, text="Print Targets", onvalue=1, offvalue=0, variable=ptargets_var)
        # self.printTargets.grid(row=0, column=5, padx=5)
        # self.cbvars["printTargets"] = ptargets_var

        # Alignment checkbox
        # pralign_var = tk.IntVar()
        # self.alignment = tk.Checkbutton(checkbox_frame, text="Alignment", onvalue=1, offvalue=0, variable=pralign_var)
        # self.alignment.grid(row=0, column=6, padx=5)
        # self.cbvars["process_aligned"] = pralign_var

        # Positional Matrix checkbox
        # pos_var = tk.IntVar()
        # self.positional_matrix = tk.Checkbutton(checkbox_frame, text="Positional Matrix", onvalue=1, offvalue=0, variable=pos_var)
        # self.positional_matrix.grid(row=0, column=7, padx=5)
        # self.cbvars["pos_matrix"] = pos_var

        # Show Positional Matrix checkbox
        # showpos_var = tk.IntVar()
        # self.show_positional_matrix = tk.Checkbutton(checkbox_frame, text="Show Positional Matrix", onvalue=1, offvalue=0, variable=showpos_var)
        # self.show_positional_matrix.grid(row=0, column=8, padx=5)
        # self.cbvars["show_pos_matrix"] = showpos_var    

        # Spacer checkbox
        spacer_var = tk.IntVar()
        self.spacer = tk.Checkbutton(checkbox_frame, text="Spacer", onvalue=1, offvalue=0, variable=spacer_var)
        self.spacer.grid(row=0, column=9, padx=5)
        self.cbvars["spacer"] = spacer_var


        
        # Update based on the checked boxes
        self.update_button = Button(self.master, text="Update", command=self.button_pressed)
        self.update_button.pack(pady=10)
        
        # Label to display the total number of sequences
        self.sequence_count_label = Label(self.master, text="")
        self.sequence_count_label.pack(pady=10)

        # Listbox to display sequence headers with a Scrollbar
        self.listbox_frame = tk.Frame(self.master)
        self.listbox_frame.pack(pady=10, padx=20, fill=tk.BOTH, expand=True)

        self.listbox_scrollbar = Scrollbar(self.listbox_frame, orient=tk.VERTICAL)
        self.listbox = Listbox(self.listbox_frame, yscrollcommand=self.listbox_scrollbar.set, height=6)
        self.listbox_scrollbar.config(command=self.listbox.yview)
        self.listbox_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.listbox.bind("<<ListboxSelect>>", self.show_sequence)


        # Create a frame for the sequence_display and its scrollbar
        self.sequence_display_frame = tk.Frame(self.master)
        self.sequence_display_frame.pack(pady=10, padx=20)

        # Create a Text widget for sequence_display
        self.sequence_display = Text(self.sequence_display_frame, height=9, width=200)
        self.sequence_display.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Create a scrollbar for the sequence_display Text widget
        self.sequence_display_scrollbar = Scrollbar(self.sequence_display_frame, orient=tk.VERTICAL, command=self.sequence_display.yview)
        self.sequence_display.config(yscrollcommand=self.sequence_display_scrollbar.set)
        self.sequence_display_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        
        # Create a frame for the tree and its scrollbar
        self.tree_frame = tk.Frame(self.master)
        self.tree_frame.pack(pady=10, padx=20, fill=tk.BOTH, expand=True)

        # Create a Treeview for the table view
        self.tree = Treeview(self.tree_frame, columns=("Name", "Description", "Length", "A", "T", "G", "C", "GC Content"), show="headings", height=7)
        self.tree.heading("Name", text="Name", command=lambda c="Name": sort_treeview(self.tree, c, False))
        self.tree.column("Name", width=300)
        self.tree.heading("Description", text="Description", command=lambda c="Description": sort_treeview(self.tree, c, False))
        self.tree.column("Description", width=600)
        self.tree.heading("Length", text="Length", command=lambda c="Length": sort_treeview(self.tree, c, False))
        self.tree.column("Length", width=50)
        self.tree.heading("A", text="A Count", command=lambda c="A": sort_treeview(self.tree, c, False))
        self.tree.column("A", width=50)
        self.tree.heading("T", text="T Count", command=lambda c="T": sort_treeview(self.tree, c, False))
        self.tree.column("T", width=50)
        self.tree.heading("G", text="G Count", command=lambda c="G": sort_treeview(self.tree, c, False))
        self.tree.column("G", width=50)
        self.tree.heading("C", text="C Count", command=lambda c="C": sort_treeview(self.tree, c, False))
        self.tree.column("C", width=50)
        self.tree.heading("GC Content", text="GC Content", command=lambda c="GC Content": sort_treeview(self.tree, c, False))
        self.tree.column("GC Content", width=50)

        # Create a scrollbar for the tree(Treeview) widget
        self.tree_scrollbar = Scrollbar(self.tree_frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.config(yscrollcommand=self.tree_scrollbar.set)

        # Pack the Treeview and its scrollbar into the frame
        self.tree_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Create a frame for the sequence_table and its scrollbar
        self.sequence_table_frame = tk.Frame(self.master)
        self.sequence_table_frame.pack(pady=10, padx=20)

        # Create a Text widget for sequence_table_display
        self.sequence_table_display = Text(self.sequence_table_frame, height=9, width=200)
        self.sequence_table_display.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Create a scrollbar for the sequence_table Text widget
        self.sequence_table_scrollbar = Scrollbar(self.sequence_table_frame, orient=tk.VERTICAL, command=self.sequence_table_display.yview)
        self.sequence_table_display.config(yscrollcommand=self.sequence_table_scrollbar.set)
        self.sequence_table_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)


        self.tree.bind("<<TreeviewSelect>>", self.show_table_sequence)
        
        # Clear button to reset selections and clear the sequence display
        self.clear_button = Button(self.master, text="Clear", command=self.clear_widget)
        self.clear_button.pack(pady=10)

        # Dictionary to store sequences from the FASTA file
        self.sequences = OrderedDict()

        # Once self.update_button has been clicked, we need to see what
        # checkboxes were checked and call the corresponding functions
        # for the selected sequence. After calling those functions, the
        # output is printed to the output table
        
    def toggle_motif_entry(self):
        # Toggle the visibility of the motif entry based on the state of the "Motif Search" checkbox
        if self.cbvars["motif"].get() == 1:
            self.motif_entry.grid()
        else:
            self.motif_entry.grid_remove()
            
    def get_header_sequence(self):
        selected_item = self.tree.focus()
        if not selected_item:
            return
        header = ">" + self.tree.item(selected_item, "values")[0] + " " + self.tree.item(selected_item, "values")[1]
        name = self.tree.item(selected_item, "values")[0]
        description = self.tree.item(selected_item, "values")[1]
        sequence = self.sequences[header]
        length = len(sequence)

        return name, description, sequence, length
    
    def dump_data(self):
        # Get the num. of pressed buttons to find the num. of columns needed
        # in the output table.
        for selected_item in self.tree.get_children():
            print(selected_item)

            header = ">" + self.tree.item(selected_item, "values")[0] + " " + self.tree.item(selected_item, "values")[1]
            name = self.tree.item(selected_item, "values")[0]
            description = self.tree.item(selected_item, "values")[1]
            sequence = self.sequences[header]
            length = len(sequence)

            print(name)

            data_processing_result = [
                (name, sequence, description)
                #(14, name, sequence, description)
            ]
            # Create a tab-delimited txt file for each table
            for table_data in data_processing_result:
                #table_name = "SEQUENCE"
                file_name = "seq.txt"
                #file_name = f"{table_name}.txt"
                #with open(file_name, "a") as file:
                with open(file_name, "a") as file:
                    file.write("\t".join(map(str, table_data))+"\n")
                print(f"File '{file_name}' created successfully")
            
    def button_pressed(self):
        # Function for when our update button is pressed.
        # Once self.update_button has been clicked, we need to see what
        # checkboxes were checked and call the corresponding functions
        # for the selected sequence. After calling those functions, the
        # output is printed to the output table
        name, description, sequence, length = self.get_header_sequence()
        self.sequence_table_display.delete(1.0, tk.END)
        codon_result = ""
        for key in self.cbvars.keys():
            if self.cbvars[key].get() == 1:
                if key == "homopolymer":
                        homopolymers = detect_homopolymer(sequence)
                        for homopolymer in homopolymers:
                            start = int(homopolymer.split("_")[0].split("-")[0])
                            end = int(homopolymer.split("_")[0].split("-")[1]) + 1
                            sequence = sequence[:start] + sequence[start:end].lower() + sequence[end:]
                            # print(start, end)
                        # print(sequence)
                        # print(homopolymers)
                        self.sequence_table_display.delete(1.0, tk.END)
                        self.sequence_table_display.insert(tk.END, sequence)
                if key == "cpg":
                    cpg_islands = cpg_island(sequence)
                    print(cpg_islands)
                    for key in cpg_islands.keys():
                        start = int(cpg_islands[key].split("_")[0].split("-")[0])
                        end = int(cpg_islands[key].split("_")[0].split("-")[1]) + 1
                        sequence = sequence[:start] + sequence[start:end].lower() + sequence[end:]
                    self.sequence_table_display.delete(1.0, tk.END)
                    self.sequence_table_display.insert(tk.END, sequence)

                if key == "motif":
                    target = self.motif_entry.get().upper()
                    if target != "":
                        motifs = motif_search(sequence, target)
                        modified_sequence = list(sequence)  # Convert the sequence to a list of characters for easy modification

                        if motifs is not None:
                            for motif_info in motifs:
                                start_end, length = motif_info.split('_')
                                start, end = map(int, start_end.split('-'))

                                for i in range(start, end + 1):
                                    modified_sequence[i] = sequence[i].lower()
                        sequence = ''.join(modified_sequence)  # Convert the modified list back to a string
                    self.sequence_table_display.delete(1.0, tk.END)
                    self.sequence_table_display.insert(tk.END, sequence)
                if key == "codon":
                    codon_dict = codon_profile(sequence)
                    codon_result = codon_profile_print(codon_dict)
                    self.sequence_table_display.delete(1.0, tk.END)
                    self.sequence_table_display.insert(tk.END, sequence)
                    self.sequence_table_display.insert(tk.END, codon_result)
                if key == "printTargets":
                    return
                if key == "process_aligned":
                    return
                if key == "pos_matrix":
                    return
                if key == "show_pos_matrix":
                    return
                if key == "spacer":
                    # print(sequence)
                    sequence = print_with_ruler(name, description, sequence, 100, True)
                    self.sequence_table_display.delete(1.0, tk.END)
                    self.sequence_table_display.insert(tk.END, sequence + f"\n\n{codon_result}")
                    # print()
                else:
                    self.sequence_table_display.delete(1.0, tk.END)
                    self.sequence_table_display.insert(tk.END, print_with_ruler(name, description, sequence, 100, False) + f"\n\n{codon_result}")

        
        # print("Button Pressed")

    def upload_file(self):
        # Open a file dialog to select the FASTA file
        fasta_file = filedialog.askopenfilename(title="Select a FASTA file",
                                                filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])
        print("================",fasta_file)

        # Display a loading message for large files
        file_size = os.path.getsize(fasta_file) / (1024 * 1024)  # File size in MB

        if file_size > 10:  # Threshold set to 10MB, can be adjusted
            loading_window = Toplevel(self.master)
            loading_label = Label(loading_window, text="Processing large file. Please wait...")
            loading_label.pack(pady=20, padx=20)
            self.master.update()
            loading_window.destroy()
        # Check if the selected file is empty
        if is_empty_file(fasta_file):
            tk.messagebox.showerror("Error", "The selected file is empty!")
            return
        if not fasta_file:
            return

        with open(fasta_file, 'r') as file:
            content = file.readlines()

        header = None
        sequence = ""
        for line in content:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    self.sequences[header] = sequence
                    self.listbox.insert(tk.END, header)
                header = line
                sequence = ""
            else:
                sequence += line
        if header:
            self.sequences[header] = sequence
            self.listbox.insert(tk.END, header)

            # Update table
            self.update_treeview()
            
            # Update the sequence count label

            loading_window.destroy()

        # Check if the FASTA file contains any sequences
        if not self.sequences:
            tk.messagebox.showerror("Error", "The selected FASTA file contains no sequences!")
            return
        self.sequence_count_label.config(text=f"Total Sequences: {len(self.sequences)}")
        

        
    def update_treeview(self):
        # Clear previous items in the Treeview
        for item in self.tree.get_children():
            self.tree.delete(item)

        # Populate the Treeview with sequence information
        for header, sequence in self.sequences.items():
            header = header[1:]
            header = header.split(" ", 1)
            name = header[0]
            description = header[1]
            length = len(sequence)
            a_count = sequence.count('A')
            t_count = sequence.count('T')
            g_count = sequence.count('G')
            c_count = sequence.count('C')
            gc_content = (g_count + c_count) / length * 100 if length > 0 else 0

            self.tree.insert("", "end", values=(name, description, length, a_count, t_count, g_count, c_count, f"{gc_content:.2f}%"))
# Modifying the ModifiedFASTAViewer class to display long sequences in a new window

class UpdatedFASTAViewer(ModifiedFASTAViewer):
    def show_sequence(self, sequence):
        selected = self.listbox.curselection()
        if not selected:
            return
        header = self.listbox.get(selected[0])
        sequence = self.sequences[header]
        
        # If sequence length exceeds 5000, open it in a new window
        if len(sequence) > 5000:
            self.show_in_new_window(sequence)
        else:
            self.sequence_display.delete(1.0, tk.END)
            self.sequence_display.insert(tk.END, sequence)

    def show_in_new_window(self, sequence):
        new_window = Toplevel(self.master)
        new_window.title("Detailed Sequence View")

        sequence_display = Text(new_window, wrap=tk.NONE, height=40, width=100)  # Increased dimensions
        sequence_scrollbar_x = Scrollbar(new_window, orient=tk.HORIZONTAL, command=sequence_display.xview)
        sequence_scrollbar_y = Scrollbar(new_window, orient=tk.VERTICAL, command=sequence_display.yview)

        sequence_display.config(xscrollcommand=sequence_scrollbar_x.set, yscrollcommand=sequence_scrollbar_y.set)

        sequence_scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X)
        sequence_scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y)
        sequence_display.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        sequence_display.insert(tk.END, sequence)

    def clear_widget(self):
        self.listbox.delete(0, tk.END)
        self.sequence_display.delete(1.0, tk.END)
        # Clear the Treeview
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.sequence_table_display.delete(1.0, tk.END)
        self.sequences.clear()
        
        # Show sequence based on selected header in table
    def show_table_sequence(self, sequence):
        selected_item = self.tree.focus()
        if not selected_item:
            return
        header = ">" + self.tree.item(selected_item, "values")[0] + " " + self.tree.item(selected_item, "values")[1]
        sequence = self.sequences[header]
        sequence = print_with_ruler(self.tree.item(selected_item, "values")[0], self.tree.item(selected_item, "values")[1], sequence, 100, False)
        # print(sequence)

        # If sequence length exceeds 5000, open it in a new window
        if len(sequence) > 5000:
            self.show_in_new_window(sequence)
        else:
            self.sequence_table_display.delete(1.0, tk.END)
            self.sequence_table_display.insert(tk.END, sequence)
            # print(self.sequence_table_display.get("1.0", "end-1c"))

    def update_table_display():
        return
    
# Function to create a MySQL database connection
def create_connection(host, user, password, database):
    connection = None
    try:
        connection = mysql.connector.connect(
            host=host,
            user=user,
            password=password,
            database=database
        )

        if connection.is_connected():
            print("Connected to MySQL Server")
            return connection

    except mysql.connector.Error as err:
        print(f"Error: {err}")
# Function to create a new database
def create_database(connection, db_name):
    try:
        cursor = connection.cursor()
        cursor.execute(f"CREATE DATABASE IF NOT EXISTS {db_name}")
        print(f"Database '{db_name}' created successfully")
    except Error as e:
        print(f"Error: {e}")

# Function to execute DDL statements for table creation
def create_tables(connection, ddl_statements):
    try:
        cursor = connection.cursor()
        for statement in ddl_statements:
            cursor.execute(statement)
        print("Tables created successfully")
    except Error as e:
        print(f"Error: {e}")

# Function to process data, extract information, and create tab-delimited text files
def process_data_and_create_txt_files(connection):
    # Replace the following with your own data processing logic
    data_processing_result = [
        ("John Doe", 25, "john.doe@example.com"),
        ("Jane Smith", 30, "jane.smith@example.com")
    ]

    # Create a tab-delimited txt file for each table
    for table_data in data_processing_result:
        table_name = "SEQUENCE"
        file_name = f"{table_name}.txt"
        with open(file_name, "w") as file:
            file.write("\t".join(map(str, table_data)))
        print(f"File '{file_name}' created successfully")

# Use this command to connect to db on ubuntu for remote access:
# mysql -u zigmonsg -p -h bio466-f15.csi.miamioh.edu


# Main script
if __name__ == "__main__":
    # Replace with your own database and user credentials
    # Remote hostname/IP address
    #host_name = "bio466-f15.csi.miamioh.edu"


    # DDL statements for table creation
    ddl_statements = [
        """
        CREATE TABLE IF NOT EXISTS SEQUENCE (
            SEQUENCE_ID INT AUTO_INCREMENT PRIMARY KEY,
            SEQUENCE_NAME VARCHAR(255),
            SEQUENCE_DESCRIPTION VARCHAR(255),
            SEQUENCE TEXT
        )
        """
        # Add more DDL statements for additional tables if needed
    ]

    # Create a connection to MySQL server
    connection = create_connection(host, user, password, database)
    #
    if connection:
        # Create the database
        create_database(connection, database)
    
        # Switch to the created database
        connection.database = database
    
        # Create tables using DDL statements
        create_tables(connection, ddl_statements)
    
        # Process data, extract information, and create tab-delimited text files
        #process_data_and_create_txt_files(connection)
    
        # Close the database connection
        connection.close()


# Import text file of sequence input into the table
# Second, you can issue the following command line described in the file seq.commandline.insert.txt
# mysqlimport -L -u zigmonsg -p bio466 seq.txt -h bio466-f15.csi.miamioh.edu
# mysqlimport -L -u zigmonsg -p bio466 seq.txt -h localhost
# mysqlimport -L -u zigmonsg -p bio466 seq.txt
# from subprocess import call
# call('mysqlimport -u zigmonsg -p bio466 seq.txt -h bio466-f15.csi.miamioh.edu', shell=True)

root = tk.Tk()
app = UpdatedFASTAViewer(root)
root.mainloop()