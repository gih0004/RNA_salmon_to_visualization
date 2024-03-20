#!/bin/bash

# Check if the number of arguments is less than three
if [ "$#" -lt 2 ]; then
  echo "Error: Insufficient arguments. Please provide two arguments, the first being refernce cdna file and the second being desirred index name."
  exit 1
fi

# Assigning the command-line arguments to variables
reference_cdna=$1  #fasta file containing transcriptome / cdna of species 
index_name=$2 #index name used for salmon


# Using the variables
echo "Reference cdna file: ${reference_cdna}"
echo "index name: ${index_name}"




#Load salmon module
module load Salmon/1.9.0-GCC-11.3.0

# make index using transcriptome fasta (used high confidence model)
salmon index -t {reference_cdna} -i {index_name}.index



for fq1 in ./.filtered_1.fq #<>This should be changed into whichever string you have last in sample names common between all samples 
do
    base="${fq1%.filtered_1.fq}"
    
#Generate transcript expression quantification files for each sample 
    salmon quant -i {index_name}.index -l A -1 ${base}.filtered_1.fq -2 ${base}.filtered_2.fq --gcBias -o ./salmon/${base}
    
done
