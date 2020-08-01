# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1 |
| 1294_S1_L008_R2_001.fastq.gz | index 1 |
| 1294_S1_L008_R3_001.fastq.gz | index 2 |
| 1294_S1_L008_R4_001.fastq.gz | read 2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. ```Your answer here```
    3. ```Your answer here```
    
## Part 2
1. Define the problem
``` a. From our original 2 FASTQ files containing our biologicaly reads (R1 and R4), we want to:
        - pair the indexes (from R2 and R3) and identify the index matches, index non-matches, and low quality/no index matches for each read
        - add the index1-index2 pair to the end of the header for each record
        - write these records into new files and count the number of records in each.
            - R1 will have:
                - 24 files for the 24 matched index pairs.
                - 1 file for unmatched/hopped index pairs.
                - 1 file for low quality index pairs.
            - R2 will have:
                - 24 files for the 24 matched index pairs.
                - 1 file for unmatched/hopped index pairs.
                - 1 file for low quality index pairs.
    b. Record 1 from each file should correspond to record 1 from each of the other files.
    c. This sequencing was run on HiSeq4000, which sequences the reverse complement of R3. Therefore, R3 needs to transformed to show this.
```
2. Describe output
```
The output files will have the index pairs added to the end of each header line, even if the indexes don't match or are of low quality.
        HEADER Index1-Index2
        SEQUENCE LINE
        +
        QUALITY SCORE LINE
```
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [4 expected output FASTQ files](../TEST-output_FASTQ).
```
See attached files.
```
4. Pseudocode
```
See attached .txt file.
```
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
