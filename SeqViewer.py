import re
import tkinter as tk
from tkinter import filedialog, Text, Listbox, Scrollbar, Toplevel, Label, Button
from collections import OrderedDict
import os
from tkinter.ttk import Treeview
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import logomaker
import pandas as pd
from NeedlemanWunsch_class_updated import Needleman_Wunsch


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
    if (
        spacer == False
    ):  # If spacer is set to False, do not add spaces between nucleotides.
        print(
            ">" + seq_name, description, "\n"
        )  # Print the sequence name and description.
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
        print(
            ">" + seq_name, description, "\n"
        )  # Print the sequence name and description.
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


class ModifiedFASTAViewer:

    def __init__(self, master):
        self.master = master
        self.master.title("FASTA Viewer")
        
        # Create a button to upload the FASTA file
        self.upload_button = Button(self.master, text="Upload FASTA", command=self.upload_file)
        self.upload_button.pack(pady=10)
        
        # List containing checkboxes
        self.cbvars = {}
        # Create a frame for checkboxes
        checkbox_frame = tk.Frame(self.master)
        checkbox_frame.pack()

        


        # Create checkboxes for settings and pack them horizontally
        spacer_var = tk.IntVar()
        self.spacer = tk.Checkbutton(checkbox_frame, text="Spacer", onvalue=1, offvalue=0, variable=spacer_var)
        self.spacer.pack(side=tk.LEFT, padx=5)
        self.cbvars["spacer"] = spacer_var

        hp_var = tk.IntVar()
        self.homopolymer = tk.Checkbutton(checkbox_frame, text="Homopolymer", onvalue=1, offvalue=0, variable= hp_var)
        self.homopolymer.pack(side=tk.LEFT, padx=5)
        self.cbvars["homopolymer"] = hp_var

        cpg_var = tk.IntVar()
        self.cpg_island = tk.Checkbutton(checkbox_frame, text="CpG Island", onvalue=1, offvalue=0, variable=cpg_var)
        self.cpg_island.pack(side=tk.LEFT, padx=5)
        self.cbvars["cpg"] = cpg_var

        # Change motif so that it takes in text input if checked
        # May just have a text field and check if it's not empty
        motif_var = tk.IntVar()
        self.motif_search = tk.Checkbutton(checkbox_frame, text="Motif Search", onvalue=1, offvalue=0, variable=motif_var)
        self.motif_search.pack(side=tk.LEFT, padx=5)
        self.cbvars["motif"] = motif_var

        codon_var = tk.IntVar()
        self.codon_profile = tk.Checkbutton(checkbox_frame, text="Codon Profile", onvalue=1, offvalue=0, variable=codon_var)
        self.codon_profile.pack(side=tk.LEFT, padx=5)
        self.cbvars["codon"] = codon_var

        prseq_var = tk.IntVar()
        self.printSeqFragment = tk.Checkbutton(checkbox_frame, text="Print Sequence Fragment", onvalue=1, offvalue=0, variable=prseq_var)
        self.printSeqFragment.pack(side=tk.LEFT, padx=5)
        self.cbvars["printSeqFragment"] = prseq_var

        ptargets_var = tk.IntVar()
        self.printTargets = tk.Checkbutton(checkbox_frame, text="Print Targets", onvalue=1, offvalue=0, variable=ptargets_var)
        self.printTargets.pack(side=tk.LEFT, padx=5)
        self.cbvars["printTargets"] = ptargets_var

        pralign_var = tk.IntVar()
        self.alignment = tk.Checkbutton(checkbox_frame, text="Aligment", onvalue=1, offvalue=0, variable=pralign_var)
        self.alignment.pack(side=tk.LEFT, padx=5)
        self.cbvars["process_aligned"] = pralign_var

        # self.process_aligned_seq = tk.Checkbutton(checkbox_frame, text="Process Aligned Sequence", onvalue=1, offvalue=0)
        # self.process_aligned_seq.pack(side=tk.LEFT, padx=5)
        
        pos_var = tk.IntVar()
        self.positional_matrix = tk.Checkbutton(checkbox_frame, text="Positional Matrix", onvalue=1, offvalue=0, variable=pos_var)
        self.positional_matrix.pack(side=tk.LEFT, padx=5)
        self.cbvars["pos_matrix"] = pos_var

        showpos_var = tk.IntVar()
        self.show_positional_matrix = tk.Checkbutton(checkbox_frame, text="Show Positional Matrix", onvalue=1, offvalue=0, variable=showpos_var)
        self.show_positional_matrix.pack(side=tk.LEFT, padx=5)
        self.cbvars["show_pos_matrix"] = showpos_var

        
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

    def button_pressed(self):
        # Function for when our update button is pressed.
        # Once self.update_button has been clicked, we need to see what
        # checkboxes were checked and call the corresponding functions
        # for the selected sequence. After calling those functions, the
        # output is printed to the output table

        # Get the num. of pressed buttons to find the num. of columns needed
        # in the output table.
        sel_functions = []
        for key in self.cbvars.keys():
            if self.cbvars[key].get() == 1:
                sel_functions.append(key)
            print(self.cbvars[key].get())

        num_functions = len(sel_functions)

        # Check the list of checkboxes to view which ones were pressed
        for func in sel_functions:
            if func == "spacer":
                # Call the specified function
                # In this case our spacer is " "
                # Still need to modify printWithRuler to output to textbox
                for header, sequence in self.sequences.items():
                    print_with_ruler(sequence, " ")
                # print()
            elif func == "homopolymer":
                # Call the specified function
                for header, sequence in self.sequences.items():
                    # Need to store this output to pass to printTargets
                    detect_homopolymer(sequence)
            elif func == "cpg":
                # Call the specified function
                for header, sequence in self.sequences.items():
                    # Need to store this output to pass to textbox
                    cpg_island(sequence)
                # print()
            elif func == "motif":
                # Call the specified function
                # Need to modify this to get text input
                print()
            elif func == "codon":
                # Call the specified function
                for header, sequence in self.sequences.items():
                    # Need to store this output to pass to codon
                    # Modify this to fit the format of this file
                    codonDict = codon_profile(sequence)
                    # print(sequence[elem])
                    print("Codon Profile:")
                    #firstLine = ">" + seq_name[elem] + " " + description[elem]
                    #print(firstLine)
                    print(sequence + "\n")

                    print("			 2nd")
                    print("         -------------------------------")
                    print("1st	   T	    C	      A      G     3rd\n")
                    # Use 3 nested loops to loop through the different codons
                    # List of nucleoties we loop through
                    nucl = ["T", "C", "A", "G"]
                    for first in range(0, len(nucl)):

                        for last in range(0, len(nucl)):
                            # String that we print out
                            if last == 0:
                                tempstr = str(first + 1) + "  	"
                            else:
                                tempstr = "        "
                            # Print out each line for each iteration of third
                            for middle in range(0, len(nucl)):
                                codon = nucl[first] + nucl[middle] + nucl[last]
                                num_matches = str(codonDict[codon])
                                tempstr += codon + "=" + " " * (3 - len(num_matches)) + num_matches + " "
                            tempstr += "   " + str(last + 1)
                            print(tempstr)
                        # Line separating the end of T codons and beginning of G codons
                        print()

                print()
            elif func == "printSeqFragment":
                # Call the specified function
                print()
            elif func == "printTargets":
                # Call the specified function
                print()
            elif func == "process_aligned":
                # Call the specified function
                print()
            elif func == "pos_matrix":
                # Call the specified function
                print()
            elif func == "show_pos_matrix":
                # Call the specified function
                print()



        print("Button Pressed")

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
        self.listbox.selection_clear(0, tk.END)
        self.sequence_display.delete(1.0, tk.END)
        
        # Show sequence based on selected header in table
    def show_table_sequence(self, sequence):
        selected_item = self.tree.focus()
        if not selected_item:
            return
        header = ">" + self.tree.item(selected_item, "values")[0] + " " + self.tree.item(selected_item, "values")[1]
        sequence = self.sequences[header]
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
            

# Display GUI
root = tk.Tk()
app = UpdatedFASTAViewer(root)
root.mainloop()