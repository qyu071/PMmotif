There are two types of data in the dataset, one is the DNA sequence dataset and the other is the motif instance information inserted in the DNA sequence dataset

Naming conventions for DNA sequence dataset files
       At-Bn-Cl-Dh-Ed-Fq-motifinf-data.txt
       The meaning of the uppercase letters in the file name
       A: The total number of DNA sequences contained in the DNA sequence dataset
       B: Base number contained in each sequence
       C: The length of the motif
       D: Start the verification layer
       E: The maximum number of mismatches allowed
       F: The number of sequences in the dataset into which motif instances are inserted
       motifinf: Motifinf for consistent sequence representation
       data: indicates that the file stores a dataset of DNA sequences

Example DNA sequence dataset file name: 3000t-200n-13l-7h-4d-600q-ATGGTAGTTATAC-data.txt

Name of the motif instance information file inserted in the DNA sequence dataset:
       At-Bn-Cl-Dh-Ed-Fq-motifinf-motif.txt

The meaning of the uppercase letters in the file name
       A: The total number of DNA sequences contained in the file
       B: Base number contained in each sequence
       C: The length of the motif
       D: Start the verification layer
       E: The maximum number of mismatches allowed
       F: The number of sequences in the dataset into which motif instances are inserted
       motifinf: Motifinf for consistent sequence representation
       motif: indicates that the file stores the information of the motif instance inserted in the DNA sequence dataset

Example motif instance information file name inserted in the DNA sequence dataset: 3000t-200n-9l-6h-2d-2100q-ACAGGGATT-motif.txt

Explanation of the meaning of the content in the document:
     
File contents
              0 181 3 CTCATACTGTT
              1 87 3 CCCTCATTGTT
              2 169 3 CCCTCACTGTA
              3 -1
              4 15 3 CCCTTACTGGG
              ...... 
              A description of what is stored in each row: Sequence Number, Insert Start Position, Hamming distance of motif and motif instance Motif instance, -1 indicates that no motif instance is inserted in the sequence.