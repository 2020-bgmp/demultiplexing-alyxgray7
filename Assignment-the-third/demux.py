#!/usr/bin/env python

### Imported methods ###
import itertools
import argparse
import gzip

### Argparse ###
def get_args():
    """
    Adds a command line option to perform this script in the terminal.
    """
    parser = argparse.ArgumentParser(description='Demultiplex a raw fastq file.')
    parser.add_argument('-r1','--read1', help='read1.fastq file location', required=True)
    parser.add_argument('-r2','--read2', help='read2.fastq file location', required=True)
    parser.add_argument('-r3','--read3', help='read3.fastq file location', required=True)
    parser.add_argument('-r4','--read4', help='read4.fastq file location', required=True)
    parser.add_argument('-cutoff','--cutoff', help='desired qscore cutoff value', required=True)
    return parser.parse_args()
args = get_args()

### Global variables ###
R1_FASTQ = args.read1
R2_FASTQ = args.read2
R3_FASTQ = args.read3
R4_FASTQ = args.read4
qscore_cutoff = int(args.cutoff)
indexes = ["GTAGCGTA","CGATCGAT","GATCAAGG","AACAGCGA","TAGCCATG","CGGTAATC","CTCTGGAT","TACCGGAT","CTAGCTCA","CACTTCAC","GCTACTCT","ACGATCAG","TATGGCAC","TGTTCCGT","GTCCTAAG","TCGACAAG","TCTTCGAC","ATCATGCG","ATCGTGGT","TCGAGAGT","TCGGATTC","GATCTTGC","AGAGTCCA","AGGATAGC"]

### Functions ###
def convert_phred(letter):
    """Converts a single character into a quality score
    single character (A) --> ASCII (65) - 33 = 32"""
    qscore = ord(letter) - 33
    return qscore

def check_qscore(qscore):
    """
    Calls on convert_phred to onverts a quality score \
        and check if the score meets the cutoff. Will \
        return a Boolean value.
    """
    for letter in qscore:
        score = convert_phred(letter)
        if score < qscore_cutoff:
            return False
    return True

def reverse_complement(sequence):
    """
    Takes a sequence string and makes a \
        reverse complement of it.
    Example: 
    sequence    reverse complement
    GTAGCGTA    TACGCTAC
    """
    complements = str.maketrans("ATGCN", "TACGN")
    revcom_sequence = sequence.translate(complements)[::-1]
    return revcom_sequence

def create_revcom_indexes():
    """
    This function will make a list of reverse complements of the indexes.
    """
    revcom_indexes = []
    for i in indexes:
        revcom_index = reverse_complement(i)
        revcom_indexes.append(revcom_index)
    return revcom_indexes
revcom_indexes = create_revcom_indexes()

def populate_matched():
    """
    This function will populate all possible index hopping combinations as a \
    tuple pair (x, y) and put them into a list.
    It will call the reverse complement function for y.
    """
    matched_dict = {}
    zipped = zip(indexes, revcom_indexes)
    matched_list = list(zipped)
    for item in matched_list:
        paired = "-".join(item)
        matched_dict.setdefault(paired,0)
        # if paired in matched_dict:
        #     matched_dict[paired] += 1
        # else:
        #     matched_dict.setdefault(paired,0)
    return matched_list, matched_dict
matched_list, matched_dict = populate_matched()
#print(matched_list)
#print(matched_dict)

def populate_unmatched():
    """
    This function will populate all possible index hopping combinations \
    as a tuple pair (x, y), where 'x' is index1 and 'y' is the reverse complement \
    of index2 and put them into a dictionary.
    """
    hopped_dict = {}
    all_indexes = [indexes, revcom_indexes]
    all_unmatched = list(itertools.product(*all_indexes))
    filtered = list(filter(lambda i: i not in matched_list, all_unmatched))
    for item in filtered:
        paired = "-".join(item)
        hopped_dict.setdefault(paired,0)
        # if paired in hopped_dict:
        #     hopped_dict[paired] += 1
        # else:
        #     hopped_dict.setdefault(paired,0)
    return hopped_dict
hopped_dict = populate_unmatched()
#print(hopped_dict)

def name_outfiles():
    """
    This function will name all the output files as strings and put \
    them into a list.
    """
    R1 = itertools.cycle(["_R1.fastq"])
    R4 = itertools.cycle(["_R4.fastq"])
    R1_files = list(zip(indexes, R1))
    pairs = ["Hopped_R1.fastq","Hopped_R4.fastq","LowQ_R1.fastq","LowQ_R4.fastq"]
    for pair1 in R1_files:
        pair1 = "".join(pair1)
        pairs.append(str(pair1))
    R4_files = list(zip(indexes, R4))
    for pair2 in R4_files:
        pair2 = "".join(pair2)
        pairs.append(str(pair2))
    return pairs
outfiles = name_outfiles()

def name_fhs():
    """
    This function will name the file handles (fh) to be used and put them in \
    to a list.
    """
    fh = itertools.cycle(["fh"])
    nums = list(range(len(outfiles)))
    y = []
    for num in nums:
        string = str(num)
        y.append(string)
        fhs = list(zip(fh, y))
        pairs = []
        for pair in fhs:
            pair = "".join(pair)
            pairs.append(pair)
    return pairs
fhs = name_fhs()

def open_files():
    """This function will open files to write for the main().
    """
    for i, item in enumerate(outfiles):
        fhs[i] = open(item,"w")
    print("Opened files")
    return

def close_files():
    """
    This function will close files written from the main().
    """
    for pos in outfiles:
        pos = (pos.split("_")[0])
        handle1 = outfiles.index(pos+"_R1.fastq")
        handle2 = outfiles.index(pos+"_R4.fastq")
        fhs[handle1].close()
        fhs[handle2].close()
    print("Closed files")
    return

### Main script ###
with gzip.open(R1_FASTQ, 'tr') as R1, gzip.open(R2_FASTQ,'tr') as R2, gzip.open(R3_FASTQ,'tr') as R3, gzip.open(R4_FASTQ,'tr') as R4:
    open_files()

    # counters
    LN = 0
    RN = 0
    lowQ_records = 0
    hopped_records = 0
    matched_records = 0

    # lists to hold each record for each file
    R1_record = []
    R2_record = []
    R3_record = []
    R4_record = []

    # read all 4 files line by line simultanously
    for R1_line, R2_line, R3_line, R4_line in zip(R1, R2, R3, R4):
        LN += 1
        R1_line = R1_line.strip()
        R2_line = R2_line.strip()
        R3_line = R3_line.strip()
        R4_line = R4_line.strip()
        
        # add lines to each record list
        if LN//4 == RN:
            R1_record.append(R1_line)
            R2_record.append(R2_line)
            R3_record.append(R3_line)
            R4_record.append(R4_line)
        else:
            R1_record.append(R1_line)
            R2_record.append(R2_line)
            R3_record.append(R3_line)
            R4_record.append(R4_line)
            RN += 1

            # index1 and index2 calculations
            index1 = R2_record[1]
            index2 = R3_record[1]
            paired = index1+"-"+index2
            rc_index2 = reverse_complement(index2)
            qscore1 = R2_record[3]
            qscore2 = R3_record[3]
            check1 = check_qscore(qscore1)
            check2 = check_qscore(qscore2)

            # add index1 and index2 onto header1 and header2
            header1 = R1_record[0] + " " + paired
            header2 = R4_record[0] + " " + paired

            # check index read quality
            if check1 and check2 == True:

                # good index read quality; check if known index read
                if index1 and rc_index2 in indexes:

                    # good index read quality; known index reads; check if paired indexes in matched_dict
                    if paired in matched_dict:
                        matched_records += 1
                        matched_dict[paired] += 1

                        # finding the index position of the index in outfiles
                        for pos in outfiles:
                            pos = (pos.split("_")[0])
                            handle1 = outfiles.index(pos+"_R1.fastq")
                            handle2 = outfiles.index(pos+"_R4.fastq")

                            # writing to the appropriate file handles for each index
                            if index1 == pos:
                                fhs[handle1].write(str(header1)+"\n"+str(R1_record[1])+"\n"+str(R1_record[2])+"\n"+str(R1_record[3])+"\n")
                                fhs[handle2].write(str(header2)+"\n"+str(R4_record[1])+"\n"+str(R4_record[2])+"\n"+str(R4_record[3])+"\n")
                                break

                    # good index read quality; known index reads; not in matched_dict; check if in hopped_dict
                    elif paired in hopped_dict:
                        hopped_records += 1
                        hopped_dict[paired] += 1
                        fhs[0].write(str(header1)+"\n"+str(R1_record[1])+"\n"+str(R1_record[2])+"\n"+str(R1_record[3])+"\n")
                        fhs[1].write(str(header2)+"\n"+str(R4_record[1])+"\n"+str(R4_record[2])+"\n"+str(R4_record[3])+"\n")
                    
                    else:
                        lowQ_records += 1
                        fhs[2].write(str(header1)+"\n"+str(R1_record[1])+"\n"+str(R1_record[2])+"\n"+str(R1_record[3])+"\n")
                        fhs[3].write(str(header2)+"\n"+str(R4_record[1])+"\n"+str(R4_record[2])+"\n"+str(R4_record[3])+"\n")

                # good index read quality; unknown index1 or rc_index2
                else:
                    lowQ_records += 1
                    fhs[2].write(str(header1)+"\n"+str(R1_record[1])+"\n"+str(R1_record[2])+"\n"+str(R1_record[3])+"\n")
                    fhs[3].write(str(header2)+"\n"+str(R4_record[1])+"\n"+str(R4_record[2])+"\n"+str(R4_record[3])+"\n")
                
            # bad index read quality
            else:
                lowQ_records += 1
                fhs[2].write(str(header1)+"\n"+str(R1_record[1])+"\n"+str(R1_record[2])+"\n"+str(R1_record[3])+"\n")
                fhs[3].write(str(header2)+"\n"+str(R4_record[1])+"\n"+str(R4_record[2])+"\n"+str(R4_record[3])+"\n")

            # clear record lists for next record
            R1_record.clear()
            R2_record.clear()
            R3_record.clear()
            R4_record.clear()

    # close output files
    close_files()

### writing outputs to output file ###
with open("outputs.txt","w") as out_file:
    # calculating percentages for index pairs
    lowQ_percent = (lowQ_records / RN) * 100
    hopped_percent = (hopped_records / RN) * 100
    matched_percent = (matched_records / RN) * 100

    # All record values table
    out_file.write("Record Type"+"\t"+"# of all"+"\t"+"Percent"+"\n")
    out_file.write("All: "+str(RN)+"\n")
    out_file.write("Matched:"+"\t"+str(matched_records)+"\t"+str(matched_percent)+"\n")
    out_file.write("Hopped:"+"\t"+str(hopped_records)+"\t"+str(hopped_percent)+"\n")
    out_file.write("Low quality/unknown:"+"\t"+str(lowQ_records)+"\t"+str(lowQ_percent)+"\n\n")
    
    # sorted table of matched indexes and counts
    out_file.write("Matched index pairs"+"\t"+"# of occurances"+"\n")
    for item in sorted(matched_dict):
        out_file.write(str(item)+"\t"+str(matched_dict[item])+"\n")
    out_file.write("\n\n")
    
    # Sorted table of hopped indexes and counts
    out_file.write("Hopped index pairs"+"\t"+"# of occurances"+"\n")
    for item in sorted(hopped_dict):
        out_file.write(str(item)+"\t"+str(hopped_dict[item])+"\n")