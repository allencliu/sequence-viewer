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



def printWithRuler(sequence, spacer):
    '''
    Print the data from the FASTA file with a spacer.
    :param sequence: The sequence we are printing here
    :param spacer: The spacer between segments of our sequence
    '''
    # Print the first line with spacers
    firstLine = "     "
    firstLineCounter = 1
    for i in range(0, 5):
        firstLine += "         " + str(firstLineCounter) + spacer
        firstLineCounter += 1
    print(firstLine)
    # Print the second line with spacers
    secondLine = "Line "
    for i in range(0, 5):
        secondLine += "1234567890" + spacer
    print(secondLine)
    # Let's use substrings to print the specific num. of nucleotides per line
    # Use a while loop to remove num. of nucleotides per line from sequence
    # We only print out 50 nucleotides per line in printWithRuler
    temp_seq = sequence
    counter = 1
    while (len(temp_seq) > 0):
        if counter < 10000 and counter >= 1000:
            curr_line = "" + str(counter) + " "
        elif counter < 1000 and counter >= 100:
            curr_line = " " + str(counter) + " "
        elif counter < 100 and counter >= 10:
            curr_line = "  " + str(counter) + " "
        else:
            curr_line = "   " + str(counter) + " "
        #curr_line = temp_seq[0, 50]
        temp_curr_line = temp_seq[0 : 50]
        while(len(temp_curr_line) > 0):
            curr_line += temp_curr_line[0 : 10] + "" + spacer
            if len(temp_curr_line) > 10:
                temp_curr_line = temp_curr_line[10:]
            else:
                temp_curr_line = ""
        print(curr_line)
        # Remove current line from temp_seq
        if len(temp_seq) > 50:
            temp_seq = temp_seq[50:]
        else:
            temp_seq = ""
        counter += 1

def detectHomopolymer(sequence):
    '''
    Function to detect the homopolymers in a DNA sequence and return
    a list with information about detected homopolymers.
    :param sequence: The DNA sequence observed here.
    :return: A list containing information about the detected
    homopolymers in the sequence (e.g., 'C 41 52 12', which represents the type,
    start position, end position and the length of detected homopolymer respectively).
    '''
    # Regex implementation, but we need to change this
    homo_list = []
    # Our regex pattern that we match
    #spattern = r'(.)\1{9,}(?:(?!\1).{0,1}\1\1{9,})*'
    pattern = r'(.)\1{9,}(?:(?!\1).{0,1}\1\1{9,})*(?:(?!\1).{0,1}\1\1{9,})?'

    # Get iteration matches using finditer()
    # Finding matches with re.finditer()
    matches = [match for match in re.finditer(pattern, sequence)]
    #matches = [match for match in re.findall(pattern, sequence)]
    #matches = re.finditer(pattern, sequence)
    #print(matches)

    # Now that we have the matches, we need to change their formatting
    # and store them in homo_list
    for match in matches:
        # Temp string to add information to
        #print(match)
        info_str = ""
        temp_span = match.span()
        match_str = match.group()
        #print(temp_span)
        info_str += str(temp_span[0]) + "_"
        info_str += str(temp_span[1]-1) + "_"
        info_str += str((temp_span[1] - temp_span[0])) + "_"
        info_str += str(match_str[0])

        #info_str = match_ch + "_" + str(start_pos) + "_" + str(end_pos) + "_" + str(len(curr_homo))
        # Print the iter matches for debugging

        #print(info_str)
        homo_list.append(info_str)

    return homo_list

def motifSearch(sequence, target):
    '''
    This function searches for matches between our DNA sequence
    and the target motif, and returns a list with information about
    those motif matches.
    :param sequence: The DNA sequence observed here.
    :param target: The target motif (a string) that we search for
    within the DNA sequence.
    :return: A list containing information about the detected motif
    matches (e.g., '41 52 12', which represents the start position,
    end position and the length of detected motif respectively).
    '''
    # List to store our mofit matches
    motif_list = []
    # Looks like we need to find EXACT matches between our target motif within
    # the sequence (no strange conditionals like for homopolymers).
    # We use a while loop sliding window technique here, and if we find a match,
    # we move our index/counter by the length of target.
    counter = 0
    while (counter < len(sequence)):
        # Check if we definitely don't have a match
        #if (counter < len(target) and sequence[counter-(len(target)-1):counter+ 1] != target):
        if (counter < len(target) and sequence[counter:counter+len(target)] != target):
            # print sequence[counter-(len(target)-1):counter+ 1] for debugging purposes
            # print(sequence[counter-(len(target)-1):counter+ 1])
            # If we don't have a match then keep sliding the window
            counter += 1
        # Check if we get a match
        elif counter == 0 and sequence[0:len(target)] == target:
            #start_pos = counter - (len(target) - 1)
            start_pos = 0
            #end_pos = counter
            end_pos = len(target) - 1
            info_str = str(start_pos) + "_" + str(end_pos) + "_" + str(len(target))
            motif_list.append(info_str)
            # Increment our counter by len(target) now that we found a match
            counter += len(target)
        elif sequence[counter-(len(target)-1):counter+ 1] == target:
            start_pos = counter-(len(target)-1)
            end_pos = counter
            info_str = str(start_pos) + "_" + str(end_pos) + "_" + str(len(target))
            motif_list.append(info_str)
            # Increment our counter by len(target) now that we found a match
            counter += len(target)
        # Else statement to increment counter normally (shouldn't be triggered/used
        # but want to avoid infinite looping
        else:
            counter += 1

    return motif_list

def printTargets(sequence, listHomo, listMotif):
    '''
    This function modifies the sequence
    :param sequence: The DNA sequence observed here.
    :param listHomo: List containing information about
    homopolymers detected in DNA sequence.
    :param listMotif: List containing information about
    motifs detected in the DNA sequence.
    :return: No return, only printing our modified sequence.
    '''
    # This function takes three inputs: a given sequence, a list of
    # detected homopolymers and a list of detected motifs. In this function,
    # you need to call your function printWithRuler() implemented in Assignment
    # No.1 to generate the sequence view with spacer and a ruler. This
    # function will not return anything out of the function.

    # Looks like we modify sequence so that all matches of homopolymers
    # and motifs in the sequence are lowercase. We do this by using the
    # indices provided in listHomo and listMotif and converting those indices
    # to lowercase.
    # Loop through listHomo and listMotif to get the indices
    for i in range(0, len(listHomo)):
        list_elems = listHomo[i].split("_")
        # (e.g., 'C 41 52 12', which represents the type, start
        # position, end position and the length of detected homopolymer respectively)
        start_pos = int(list_elems[0])
        end_pos = int(list_elems[1])
        # Now convert elements within this subsequence to lowercase
        # print(start_pos)
        # print(end_pos)
        # print(sequence)
        # print(sequence[start_pos:end_pos+1])
        lower_sub = sequence[start_pos:end_pos + 1].lower()
        # Split up the sequence by strings before the subsequence and after
        before = sequence[:start_pos]
        after = sequence[end_pos + 1:]
        # Add the lower_sub in between the before and after
        sequence = before + lower_sub + after

    # Loop through listMotif to get the indices
    for i in range(0, len(listMotif)):
        list_elems = listMotif[i].split("_")
        # (e.g., '41 52 12',  which represents the start position,
        # end position and the length of detected motif respectively)
        start_pos = int(list_elems[0])
        end_pos = int(list_elems[1])
        # Now convert elements within this subsequence to lowercase
        #print(start_pos)
        #print(end_pos)
        #print(sequence)
        #print(sequence[start_pos:end_pos+1])
        lower_sub = sequence[start_pos:end_pos + 1].lower()
        # Split up the sequence by strings before the subsequence and after
        before = sequence[:start_pos]
        after = sequence[end_pos + 1:]
        # Add the lower_sub in between the before and after
        sequence = before + lower_sub + after

    # Once we converted everything to lowercase, now we print using printWithRuler()
    # (not sure about the spacer parameter, but I can just default to " ")
    printWithRuler(sequence, "")

def CpGIsland(sequence):
    '''
    Fina all instances of
    :param sequence: The DNA sequence we are given
    :return: A dictionary with index locations of CpGIslands
    as keys and the specific CpGIsland as values
    '''
    # We need to use for loops (not regex) to
    # obtain a dictionary of keys of detected CpGIsland
    # numbers
    # Print this for debugging purposes
    #print("Coming inside ", sequence)
    # Dictionary we return of CpGIslands as
    # values and index locations as keys
    cpgDict = {}
    # Counter for # of successive CpG pairs
    counter = 0
    # Loop through all elements in steps of 1
    for i in range(0, len(sequence)-1, 1):
        if (sequence[i] + sequence[i+1]) == "CG":
            counter += 1
        # We need to add to our dict if we get 6 pairs
        elif counter >= 12:
            key = str(i-(counter+1)) + "-" + i-1 + "_" + counter
            cpgDict[key] = sequence[i-(counter+1), i]
        else:
            counter = 0

    #print("End of CpGIsland")
    return cpgDict

def translation(sequence):
    '''
    Takes DNA sequence as an input and returns a protein
    sequence (string) after translation.
    :param sequence: The DNA sequence we are given
    :return: A protein sequence (string)
    '''
    # Use the translation function provided in first class
    # of this course, but make sure to use nucleotide “T”
    # to replace “U” because of DNA sequence mentioned here.

    # File name: translation.py
    trans_dic = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
                 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
                 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
                 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
                 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
                 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                 'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
                 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                 'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
                 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

    # Not sure if I need to make additional changes
    # here to reflect translation. I don't think
    # I need to reverse the sequence or change out
    # base pairs with their complements here
    # (happens in transciption?).
    dnaseq = sequence
    rnaseq = dnaseq.replace('T', 'U')
    # Replace N base with empty space
    rnaseq = rnaseq.replace('N', '')
    protein = []
    for i in range(0, len(rnaseq)-2, 3):
        protein.append(trans_dic[rnaseq[i:i+3]])
    #print("".join(protein))
    # We need to replace U with T since we
    # are returning a DNA sequence and not RNA
    #result_dna = "".join(protein)
    #result_dna = result_dna.replace('U', 'T')
    return "".join(protein)

def complement(seq_fragment):
    '''
    Helper function to get the complement of
    a sequence fragment.
    :param seq_fragment: The sequence fragment in question.
    :return: The complement of the sequence fragment.
    '''
    comp = seq_fragment.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c')
    #revcomp = rev.upper()[::-1]
    comp = comp.upper()
    return comp

def translation6Frame(sequence):
    '''
    :param sequence: The DNA sequence we are given
    :return: A tuple of protein sequences produced by
    the sequence.
    '''
    # Next, you need to write a new function
    # translation6Frame(), which takes DNA sequence
    # as the input and returns 6 protein sequences
    # as a tuple. Of course, within translation6Frames(),
    # you need to call translation6Frame()
    # (or do I call translation() 6 times?)
    # six times to get these 6 protein sequences.

    # Split the DNA sequence into thirds. We call translation
    # on each of these thirds and add the result to listProtein.
    # Reverse the DNA sequence, and repeat the process on the
    # thirds of the reversed sequence. Then we convert listProtein
    # to a tuple and return it

    listProtein = []
    # Split the sequence into 3 substrings
    #third = int(len(sequence)/3)
    #listProtein.append(translation(sequence[0:third]))
    #listProtein.append(translation(sequence[third:third*2]))
    #listProtein.append(translation(sequence[third*2:len(sequence)]))

    # Reverse the sequence
    rev_seq = sequence[::-1]
    #rev_seq = complement(rev_seq)
    #listProtein.append(translation(rev_seq[0:third]))
    #listProtein.append(translation(rev_seq[third:third * 2]))
    #listProtein.append(translation(rev_seq[third*2:len(sequence)]))
    listProtein.append(translation(rev_seq))
    listProtein.append(translation(rev_seq[1:]))
    listProtein.append(translation(rev_seq[2:]))

    listProtein.append(translation(sequence))
    listProtein.append(translation(sequence[1:]))
    listProtein.append(translation(sequence[2:]))
    # Question: What part of the sequence are we pulling
    # the proteins from, especially since are only returning
    # a tuple of 6 protein sequences? Ask Liang about this

    # Convert listProtein to a tuple
    listProtein = tuple(listProtein)
    return listProtein


# Need to revise printSeqFragment, one of the fragments on the top
# or bottom is incorrect
def printSeqFragment(seq_fragment, start, end):
    '''

    :param seq_fragment: The sequence fragment to
    be printed.
    :param start: The starting position we print from.
    :param end: The ending position we stop printing
    the sequence fragment from.
    :return:
    '''
    # What does the output (letters for proteins) mean
    # for printSeqFragment? How are they supposed to be
    # structured? Ask Dr. Liang about this.
    # It seems like printSeqFragment calls translation6Frame

    # Sounds like we want to actually use the gene complement
    # on the top and bottom rows instead of reversing it

    # Look at reverse of sequence fragment instead?
    comp = complement(seq_fragment[start:end])
    # Tuple storing the 6 protein frames
    proteins = translation6Frame(seq_fragment[start:end])
    #proteins = translation6Frame(comp)

    # The complement sequence and protein frames are at the top
    # Each of the 3 protein layers on top is each of the frames
    # gathered from translation6frame corresponding to the complement/
    # original sequence passed (first 3 frames in the tuple).

    # Print the first 3 protein fragments
    # May need to reformat how these are printed
    # print(' '.join(intersection_lp))
    #print(proteins[0])
    #print(proteins[1])
    #print(proteins[2])
    print("  " + "  ".join(proteins[0]))
    print("   " + "  ".join(proteins[1]))
    print(" " + "  ".join(proteins[2]))
    # The selected fragment has a length of 42 nucleotides:
    # Now we print the first and last 14 chars of the sequence
    #first_pos = str(seq_fragment[elem][start: start + 14]) + "..." \
                #+ str(seq_fragment[elem][end - 15: end + 1])
    rev_fragment = seq_fragment[::-1]
    # Print out complementary sequence instead?
    first_pos = str(comp)
    #first_pos = str(seq_fragment[start:end])
    print(first_pos)
    second_pos = "|" * len(first_pos)
    third_pos = "<" + str(start) + "-" * (len(first_pos) - (len(str(start)) + len(str(end)) + 2)) \
                + str(end) + ">"
    print(second_pos)
    print(third_pos)
    # AGGATATAGATAAGATAATAGACCACCCATAGA
    # |||||||||||||||||||||||||||||||||
    # <58---------------------------99>


    # The reversed sequence is at the bottom and the protein frames
    # are below it. Each of the 3 protein layers on bottom is each of the frames
    # gathered from translation6frame corresponding to the reversed sequence
    # passed (last 3 frames in the tuple).



    # The selected fragment has a length of 42 nucleotides:
    # Now we print the first and last 14 chars of the sequence
    bottom_first_pos = str(seq_fragment[start:end])
    #bottom_first_pos = str(rev_fragment[start:end])
    print(bottom_first_pos)
    bottom_second_pos = "|" * len(first_pos)
    print(bottom_second_pos)
    # |||||||||||||||||||||||||||||||||
    # AGGATATAGATAAGATAATAGACCACCCATAGA

    # Print the last 3 protein fragments
    # May need to reformat how these are printed
    #print(proteins[3])
    #print(proteins[4])
    #print(proteins[5])
    print(" " +"  ".join(proteins[3]))
    print("  " + "  ".join(proteins[4]))
    print("   " + "  ".join(proteins[5]))

def alignment(seq_name_list, seq_list, score_list):
    '''
    A function to perform sequence alignment on multiple DNA
    sequences using the Needleman Wunsch algorithm.
    :param seq_name_list: A list of DNA sequence names.
    :param seq_list: A list of DNA sequences.
    :param score_list: A list of the scoring system we use for
    match, mismatch, and gap, respectively.
    :return: Does not return, only prints output (by calling the
    process_aligned_seq() function).
    '''
    # For sequence alignment, you need to implement a new function
    # alignment(seq_name_list, seq_list), which takes a list of
    # sequence names and a list of sequences as the input. This
    # function will not return anything. In this function, you will
    # use NeedlemanWunsch_class module to conduct pairwise alignments multiple
    # times. You will do the sequence alignment between the first
    # sequence and the second, the first and the third, and the first and fourth, etc.
    # For two aligned sequences that you get for each pair-wise alignment,
    # you will call a function named process_aligned_seq(aligned_seq1, aligned_seq2)
    # several times to print the alignments and other information including CIGAR string.

    # This function performs sequence alignment using Needleman Wunsch algorithm
    # and alignments and other relevant info are printed in process_aligned_seq
    # We will do pairwise alignment using 2 for loops
    #for i in range(0, len(seq_list)):
    # Concretely, you use the first selected sequence as the
    # reference, and do the pairwise alignment with the second,
    # third, fourth one, until the last one in the selected
    # sequences if more than 4 are selected
    # In this case, we just set i = 0 since we compare first sequence
    # to all the others
    i = 0
    for j in range(0, len(seq_list)):
        # Perform pairwise alignment for different sequences
        if i != j:
            # Perform Needleman Wunsch on our sequences
            # May need to change scores for match, mismatch, and gap
            #objA = Needleman_Wunsch(seq_list[i], seq_list[j], 2, -1, -2)
            #objA = Needleman_Wunsch(seq_list[i], seq_list[j], score_list[0],
            #                        score_list[1], score_list[2])
            objA = Needleman_Wunsch(seq_list[i], seq_list[j], score_list[0],
                                    score_list[1], score_list[2])
            aligned_seq1, aligned_seq2 = objA.give_final_result()

            # Call process_aligned_seq() with our aligned sequences
            # Should also pass in sequence names here
            process_aligned_seq(aligned_seq1, aligned_seq2, seq_name_list[i], seq_name_list[j])


def process_aligned_seq(aligned_seq1, aligned_seq2, name_seq1, name_seq2):
    '''

    :param aligned_seq1:
    :param aligned_seq2:
    :return:
    '''
    # Within this function, you will generate the middle line of
    # the alignment and CIGAR string.  You also need to count the
    # matches, mismatches, insertions or deletions and compute their ratios.
    # All the print statements for generating the alignment-related printout
    # should be implemented within this function.

    # Output should look like this:
    # You need to display this translation for all selected 4 sequences.
    # [2] Pairwise sequence alignments for the sequence fragments of the selected sequences:
    # ENSEMBL008   AGGATATAGATAAGATAGACCACCCATAGA
    #              |||| |||||||||||||| |||||:|| |
    # ENSEMBL009   AGGA-ATAGATAAGATAGA-CAGGCGTATA
    #
    # Identities=26/30 (86.7%) Similarity=27/30 (90.0%) Gaps=2/30 (6.6%) CIGAR = 4M1D14M1D10M
    # #  You need to display this pairwise sequence alignments for all selected sequence using the
    # #  first one as the reference.
    # Check for purines or pyrimidines for the : part of alignment
    # R: purine (A or G)       Y: pyrimidine (C or T)

    # Count Sequence Identity (number of matches/identical nucleotides)
    seq_identity = 0
    # Count Sequence Similarity (number of matches and similar nucleotides like
    # similar purines, pyrimidines)
    seq_similarity = 0
    # Count number of gaps
    num_gaps = 0
    # Loop through both sequences and build a string of matches, mismatches, gaps, and :
    # Assumes both aligned sequences are same length
    al_str = ""
    for i in range(0, len(aligned_seq1)):
        # Print for debugging
        #print(aligned_seq1[i], aligned_seq2[i])

        # Check for matches
        if (aligned_seq1[i] == aligned_seq2[i]):
            al_str += "|"
            seq_identity += 1
            seq_similarity += 1
        # Check for gaps
        elif (aligned_seq1[i] == "-" or aligned_seq2[i] == "-"):
            al_str += " "
            num_gaps += 1
        # Check for purines or pyrimidines matches (but overall mismatches)
        elif (aligned_seq1[i] == "A" and aligned_seq2[i] == "G") or \
            (aligned_seq1[i] == "G" and aligned_seq2[i] == "A") or \
               (aligned_seq1[i] == "C" and aligned_seq2[i] == "T") or \
            (aligned_seq1[i] == "T" and aligned_seq2[i] == "C"):
            al_str += ":"
            seq_similarity += 1
        # Otherwise, it is a mismatch
        else:
            al_str += " "


        cigar_string = ""
        start = 0
        end = 0

        # Iterate over both sequences to identify matches and gaps
        for base1, base2 in zip(aligned_seq1, aligned_seq2):
            if base1 == base2:
                end += 1
            else:
                # If there's a mismatch, or a gap (insertion or deletion), append to the CIGAR string
                if start != end:
                    # In this case we have a match
                    cigar_string += str(end - start) + "M"
                start = end + 1
                end = start
                if base1 == '-':
                    # In this case we have insertion in sequence1
                    cigar_string += "1I"
                elif base2 == '-':
                    # In this case we have deletion in sequence1
                    cigar_string += "1D"

        # Add any remaining matches at the end of the sequences
        if start != end:
            cigar_string += str(end - start) + "M"
        cigar = cigar_string

    # Now we print everything out
    first_line = name_seq1 + "   " + aligned_seq1
    second_line = len(name_seq1) * " " + "   " + al_str
    third_line = name_seq2 + "   " + aligned_seq2
    print()
    print(first_line)
    print(second_line)
    print(third_line)
    # ENSEMBL008   AGGATATAGATAAGATAGACCACCCATAGA
    #              |||| |||||||||||||| |||||:|| |
    # ENSEMBL009   AGGA-ATAGATAAGATAGA-CAGGCGTATA
    print()
    # Now print out additional info
    # Get CIGAR format
    ident = (seq_identity/len(aligned_seq1)) * 100
    ident = round(ident, 2)
    similarity = (seq_similarity/len(aligned_seq1)) * 100
    similarity = round(similarity, 2)
    gaps = (num_gaps/len(aligned_seq1)) * 100
    gaps = round(gaps, 2)
    fifth_line = "Identities=" + str(seq_identity) + "/" + str(len(aligned_seq1)) + " (" + \
                 str(ident) + "%)" + " Similarity=" + str(seq_similarity) + \
                 "/" + str(len(aligned_seq1)) + " (" + str(similarity) + \
                 "%)" + " Gaps=" + str(num_gaps) + "/" + str(len(aligned_seq1)) + " (" + \
                 str(gaps) + "%)" + " CIGAR = " + cigar
    print(fifth_line)
    # Identities=26/30 (86.7%) Similarity=27/30 (90.0%) Gaps=2/30 (6.6%) CIGAR = 4M1D14M1D10M


def positional_matix(seq_list):
    '''
    Function to generate a positional matrix of nucleotides
    from a list of DNA sequences.
    :param seq_list: List of DNA sequences
    :return: A 2-dimensional NumPy array that represents a
    positional matrix of nucleotides.
    '''
    # The function positional_matix(seq_list) will take a list of
    # sequences as the only input. It will return a 2-dimensional
    # numpy array that represents the positional matrix. In the lecture on
    # Numpy array, we provide a good code example from which you can learn from.
    # Dictionary to store indices for nucleotides
    nuc_dict = {"A": 0, "T": 1, "G": 2, "C": 3}
    # Initialize our positional matrix
    pos_matrix = np.zeros((4, len(seq_list[0])), dtype=int)
    # Loop through DNA sequences in DNA list
    for i in range(0, len(seq_list)):
        # Loop through bases within the current DNA sequence
        # to add to positional matrix. We use enumerate() to get
        # the index and element
        for inx, nucl in enumerate(seq_list[i]):
            pos_matrix[nuc_dict[nucl]][inx] += 1

    return pos_matrix

def show_positional_matix(my_array):
    '''
    Function to print the contents of the positional nucleotide matrix
    and save a stacked bar chart of the nucleotides.
    :param my_array: The 2-d NumPy array for positional nucleotide
    matrix generated in positional_matrix
    :return: Does not return anything. Instead prints the positional
    nucleotide matrix and generates and saves an over-stacked bar chart.
    '''
    # .  It will not return anything. In this function,
    # you need to print the positional nucleotide matrix as shown
    # in Screen 4. You also need to generate an over-stacked bar
    # chart and save it as an image file named “PNM.jpg”.
    # nuc_dict = {"A": 0, "T": 1, "G": 2, "C": 3}
    nuc_list = ["A", "T", "G", "C"]
    first_line = "Position   "
    # Loop through all positions in positional matrix
    for i in range(0, len(my_array[0])):
        if i < 10:
            first_line += str(i+1) + "   "
        elif i < 100:
            first_line += str(i + 1) + "  "

    second_line = "--------   "
    for i in range(0,len(my_array[0])):
        if i < 10:
            second_line += "--  "
        elif i < 100:
            second_line += "--  "

    print(first_line)
    print(second_line)
    # Loop through each row in positional matrix
    for inx, row in enumerate(my_array):
        curr_line = ""
        # Loop through each column in positional matrix
        for i in range(0, len(row)):
            if i == 0:
                curr_line += str(nuc_list[inx]) + " " * 10
                curr_line += str(row[i]) + "   "
            # elif i < 10:
            #     curr_line += str(row[i]) + "   "
            # elif i < 100:
            #     curr_line += str(row[i]) + "  "
            # Assume we have less than 1000 nucleotides in our sequence
            else:
                curr_line += str(row[i]) + "   "
        print(curr_line)

    # Position	1	2	3	4	5	6	7	8	9   ……
    # --------	--	--	--	--	--	--	--	--	--
    # A	4	0	4	0	0	2	0	0	2
    # T	0	1	0	3	0	0	4	0	0
    # G	0	3	0	1	4	0	0	0	0
    # C	0	0	0	0	0	2	0	4	2

    # Generate plots using matplotlib
    # This stacked bar chart uses position for the x-axis (colum #)
    # and percentage of occurrences as the y-axis.
    plt.rcParams["figure.figsize"] = (10, 10)
    index = np.arange(0, len(my_array[0]))
    width = 0.40

    # Change my_array to percentage, divide by 4 (total num of nucleotides
    # in each column), and then multiply by 100
    # my_array = np.divide(my_array, 4)
    # my_array = my_array * 100
    #print(my_array)
    # Get nucleotide counts for each position
    a_counts = my_array[0]
    t_counts = my_array[1]
    g_counts = my_array[2]
    c_counts = my_array[3]

    plt.bar(index, a_counts, width, color='blue', label='A')
    #plt.bar(index, t_counts, width, color='orange', label='Girls Grades', bottom=a_counts)
    plt.bar(index, t_counts, width, color='orange', label='T', bottom=a_counts)
    plt.bar(index, g_counts, width, color='grey', label='G', bottom=a_counts + t_counts)
    plt.bar(index, c_counts, width, color='yellow', label='C', bottom=a_counts + t_counts + g_counts)
    plt.title("Positional Nucleotide Matrix")
    #plt.ylabel("Grades")
    #plt.xlabel("Schools")
    plt.xticks(index)
    plt.legend(loc='best')
    #plt.show()
    plt.savefig("PNM.jpg")
    plt.show()

    # We also add the bonus plot here. In order to get this bonus points,
    # you need to generate a logo graph using
    # https://github.com/jbkinney/logomaker/blob/master/logomaker
    # or any other python module that you can find.
    # This appears to be the stacked bar chart from before, but
    # the stacked bars are replaced with stacked letters of nucleotides

    import logomaker
    # Make my_array into a dataframe to be used with logomaker
    pos_df = pd.DataFrame(my_array.T)
    pos_df.columns = ["A", "T", "G", "C"]
    #print(pos_df)
    # print(logo)
    logo = logomaker.Logo(pos_df, font_name='Arial Rounded MT Bold')
    logo.fig.savefig("Logo.jpg")
    logo.fig.show()

def codonProfile(sequence):
    '''
    Function to create a dictionary of 64 codons as keys
    and their frequencies as values.
    :param sequence: The DNA sequence we are given.
    :return:
    '''
    # Initialize the dictionary with 64 codons as the keys
    # with the values of 0 first. We achieve this using a loop.
    codonDict = {}
    # I think I start this with initializing a list of all the
    # codons, loop through this list to make them keys in our
    # dictionary, setting those values to 0.
    codonList = []

    nucl = ["A", "T", "G", "C"]
    # Using nested for loops to generate all codon combos
    for first in range(0, len(nucl)):
        for second in range(0, len(nucl)):
            for third in range(0, len(nucl)):
                codonList.append(str(nucl[first] + nucl[second]+ nucl[third]))

    for codon in codonList:
        codonDict[codon] = 0
    # Now we need to fill codonDict with the proper frequencies
    # of codons within our sequence
    # Loop by 3 codons at a time
    #print(sequence)
    # Don't jump by 3, we want to capture all frames
    for i in range(0, len(sequence) - 2, 1):
        currCodon = str(sequence[i]) + str(sequence[i+1]) + str(sequence[i+2])
        #print(currCodon)
        if currCodon in codonDict.keys():
            codonDict[currCodon] += 1
        if (i + 3) > len(sequence) - 1:
            break

    return codonDict


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
        self.spacer = tk.Checkbutton(checkbox_frame, text="Spacer", onvalue=1, offvalue=0)
        self.spacer.pack(side=tk.LEFT, padx=5)

        self.homopolymer = tk.Checkbutton(checkbox_frame, text="Homopolymer", onvalue=1, offvalue=0)
        self.homopolymer.pack(side=tk.LEFT, padx=5)

        self.cpg_island = tk.Checkbutton(checkbox_frame, text="CpG Island", onvalue=1, offvalue=0)
        self.cpg_island.pack(side=tk.LEFT, padx=5)

        self.motif_search = tk.Checkbutton(checkbox_frame, text="Motif Search", onvalue=1, offvalue=0)
        self.motif_search.pack(side=tk.LEFT, padx=5)

        self.codon_profile = tk.Checkbutton(checkbox_frame, text="Codon Profile", onvalue=1, offvalue=0)
        self.codon_profile.pack(side=tk.LEFT, padx=5)

        self.printSeqFragment = tk.Checkbutton(checkbox_frame, text="Print Sequence Fragment", onvalue=1, offvalue=0)
        self.printSeqFragment.pack(side=tk.LEFT, padx=5)

        self.printTargets = tk.Checkbutton(checkbox_frame, text="Print Targets", onvalue=1, offvalue=0)
        self.printTargets.pack(side=tk.LEFT, padx=5)

        self.alignment = tk.Checkbutton(checkbox_frame, text="Aligment", onvalue=1, offvalue=0)
        self.alignment.pack(side=tk.LEFT, padx=5)

        self.process_aligned_seq = tk.Checkbutton(checkbox_frame, text="Process Aligned Sequence", onvalue=1, offvalue=0)
        self.process_aligned_seq.pack(side=tk.LEFT, padx=5)

        self.positional_matrix = tk.Checkbutton(checkbox_frame, text="Positional Matrix", onvalue=1, offvalue=0)
        self.positional_matrix.pack(side=tk.LEFT, padx=5)

        self.show_positional_matrix = tk.Checkbutton(checkbox_frame, text="Show Positional Matrix", onvalue=1, offvalue=0)
        self.show_positional_matrix.pack(side=tk.LEFT, padx=5)
        
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
            if self.cbvars[key] == True:
                sel_functions.append[key]

        num_functions = len(sel_functions)

        # Check the list of checkboxes to view which ones were pressed
        for func in sel_functions:
            if func == "spacer":
                # Call the specified function
                # In this case our spacer is " "
                # Still need to modify printWithRuler to output to textbox
                for header, sequence in self.sequences.items():
                    printWithRuler(sequence, " ")
                # print()
            elif func == "homopolymer":
                # Call the specified function
                for header, sequence in self.sequences.items():
                    # Need to store this output to pass to printTargets
                    detectHomopolymer(sequence)
            elif func == "cpg":
                # Call the specified function
                for header, sequence in self.sequences.items():
                    # Need to store this output to pass to textbox
                    CpGIsland(sequence)
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
                    codonDict = codonProfile(sequence)
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