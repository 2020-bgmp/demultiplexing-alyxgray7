#!/usr/bin/env python

import argparse
def get_args():
    """
    Add a command line option to perform this script in the terminal.
    """
    parser = argparse.ArgumentParser(description='Create a desired length for a list array. Specify the desired file.')
    parser.add_argument('-f','--file', help='Desired file', required=True)
    parser.add_argument('-l','--length', type=int, help='Length of sequence read')
    parser.add_argument('-o','--output', help='Read identification (R1, R2, R3, or R4)', required=True)
    return parser.parse_args()
args = get_args()

file_in = args.file
seq_length = args.length
read = args.output

# FUNCTION - CONVERT TO QUALITY SCORE
def convert_phred(letter):
    """Converts a single character into a quality score
    single character (A) --> ASCII (65) - 33 = 32"""
    qscore = ord(letter) - 33
    return qscore

# FUNCTION - INITIALIZE A LIST
def init_list(array, value=0.0):
    '''This function takes an empty list and will populate it with
    the value passed in "value".'''
    for i in range(seq_length):
        array.append(value)
    return array

# initiate variables
score_sums = []
score_sums = init_list(score_sums)

# FUNCTION - POPULATE LIST WITH QUALITY SCORE SUMS
def populate_sums(file_in):
    """This function will sum the quality scores at each base position."""
    
    # variables
    LN = 0
    
    import gzip
    with gzip.open(file_in,"tr") as fh1:
        for line in fh1:
            LN += 1
            
            # progress print statement
            if LN%100000000 == 0:
                print("Finished: ", LN)

            # extract the quality score lines
            if LN%4 == 0:
                letter_counter = 0
                
                # converting lines to score sums
                for letter in line.rstrip("\n"):
                    score = convert_phred(letter)
                    score_sums[letter_counter] += score
                    letter_counter += 1
    return score_sums, LN  

score_sums, LN = populate_sums(file_in)

# FUNCTION - CALCULATE QSCORE MEANS
def calculate_mean(file_in):
    """This function will calculate the mean quality score at each base position."""

    # variables
    mean_scores = []

    # divide the sum at each position by the number of reads
    for running_sum in score_sums:
        records = LN/4
        average = running_sum/records
        mean_scores.append(average)
    return mean_scores

mean_scores = calculate_mean(file_in)

# FIGURE MAKING
import matplotlib.pyplot as plt
x = range(seq_length)
y = mean_scores

plt.bar(x,y)
plt.title("Mean Quality Score Distribution_{}".format(read))
plt.xlabel('Base Position')
plt.ylabel('Mean Quality Score')
plt.savefig("Distribution_{}".format(read))